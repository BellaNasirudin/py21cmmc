"""
A module providing Core Modules for cosmoHammer. This is the basis of the plugin system for py21cmmc.
"""
import os
import numpy as np
from . import wrapper as p21c
import string


class Core21cmFastModule:
    storage_options = {
        "DATADIR": "MCMCData",
        "KEEP_ALL_DATA": True,
        "KEEP_GLOBAL_DATA": True
    }

    # The following disables the progress bar for perturbing and ionizing fields.
    disable_progress_bar = True

    def __init__(self, parameter_names,
                 storage_options={},
                 box_dim={}, flag_options={}, astro_params={}, cosmo_params={}):

        # SETUP variables
        self.parameter_names = parameter_names

        self.storage_options.update(**storage_options)

        # Save the following only as dictionaries (not full Structure objects), so that we can avoid
        # weird pickling errors.
        self._box_dim = box_dim
        self._flag_options = flag_options
        self._astro_params = astro_params
        self._cosmo_params = cosmo_params

    def setup(self):
        print("Using Core21cmFAST Module")

        self._regen_init = False
        self._write_init = True
        self._modifying_cosmo = False
        # The following ensures that if we are changing cosmology in the MCMC, then we re-do the init
        # and perturb_field parts on each iteration.
        for p in self.parameter_names:
            if p in p21c.CosmoParamStruct._defaults_.keys():
                self._flag_options['GenerateNewICs'] = True
                self._write_init = False
                self._regen_init = True
                self._modifying_cosmo = True

        # Here we create the init boxes and perturb boxes, written to file.
        # If modifying cosmo, we don't want to do this, because we'll create them
        # on the fly on every iteration. We don't need to save any values in memory because
        # they will be read from file.
        if not self._modifying_cosmo:
            p21c.run_21cmfast(
                self._flag_options['redshifts'],
                self._box_dim,
                self._flag_options,
                run_ionize=False,
                write=True, regenerate=self._regen_init,
                progress_bar=not self.disable_progress_bar
            )

    def __call__(self, ctx):
        AstroParams, CosmoParams = self._update_params(ctx)
        output = self._run_21cmfast(AstroParams, CosmoParams)

        ctx.add('output', output)
        ctx.add('box_dim', self.box_dim())
        ctx.add('cosmo_params', self.cosmo_params(**CosmoParams))
        ctx.add("astro_params", self.astro_params(**AstroParams))
        ctx.add("flag_options", self.flag_options())

    def box_dim(self, **kwargs):
        dct = dict(self._box_dim, **kwargs)
        return p21c.BoxDimStruct(**dct)

    def cosmo_params(self, **kwargs):
        dct = dict(self._cosmo_params, **kwargs)
        return p21c.CosmoParamStruct(**dct)

    def astro_params(self, **kwargs):
        dct = dict(self._astro_params, **kwargs)
        return p21c.AstroParamStruct(self._flag_options.get("INHOMO_RECO", False), **dct)

    def flag_options(self, **kwargs):
        dct = dict(self._flag_options, **kwargs)
        # TODO: it is bad that we call astr_params() here without kwargs...
        return p21c.FlagOptionStruct(self.astro_params().Z_HEAT_MAX, self.astro_params().ZPRIME_STEP_FACTOR,
                                     **dct)

    def _update_params(self, ctx):
        """
        Update all the parameter structures which get passed to the driver, for this iteration.

        Parameters
        ----------
        ctx : context
        """
        params = ctx.getParams()

        # Copy the main parameter structures (can't write to them otherwise each walker will over-write each other)
        #        AstroParams = deepcopy(self._astro_params)
        #        CosmoParams = deepcopy(self._cosmo_params)

        apkeys = self.astro_params().c_initializer
        cpkeys = self.cosmo_params().c_initializer

        AstroParams = {}
        CosmoParams = {}

        # Update the Astrophysical/Cosmological Parameters for this iteration
        for k in params.keys:
            if k in apkeys:
                AstroParams[k] = getattr(params, k)
            elif k in CosmoParams:
                CosmoParams[k] = getattr(params, k)
            else:
                raise ValueError("Key %s is not in AstroParams or CosmoParams " % k)

        if self._regen_init:
            # A random number between 1 and 10^12 should be sufficient to randomise the ICs
            CosmoParams['RANDOM_SEED'] = int(np.random.uniform(low=1, high=1e12, size=1))

        return AstroParams, CosmoParams

    def _run_21cmfast(self, AstroParams, CosmoParams):
        """
        Actually run the 21cmFAST code.

        Parameters
        ----------

        Returns
        -------
        lightcone : Lightcone object.
        """
        if self._modifying_cosmo:
            return p21c.run_21cmfast(
                self._flag_options['redshifts'],
                self._box_dim,
                self._flag_options,
                AstroParams, CosmoParams,
                self._write_init, self._regen_init
            )[0]
        else:
            return p21c.run_21cmfast(
                self._flag_options['redshifts'],
                self._box_dim,
                self._flag_options,
                AstroParams, CosmoParams,
                self._write_init, self._regen_init,
            )[0]

    def _store_data(self, lightcone, random_ids):
        """
        Write desired data to file for permanent storage.

        Currently unused, but ported from Brad's verison.

        Parameters
        ----------
        lightcone : the lightcone object

        """

        StringArgument_other = "%s_%s" % random_ids
        StoredStatisticalData = []
        StoredFileLayout = []

        separator_column = "\t"

        if not self.storage_options['KEEP_GLOBAL_DATA']:
            if not self.use_lightcone:

                for i in range(len(self.redshift)):

                    if self.storage_options['KEEP_ALL_DATA']:

                        if i == 0:
                            StoredStatisticalData.append(lightcone.k)
                            StoredFileLayout.append("{%i}" % (i))

                        StoredStatisticalData.append(lightcone.power_spectrum[i])
                        StoredFileLayout.append("{%i}" % (i + 1))

            if self.storage_options['KEEP_ALL_DATA']:

                StoredFileLayout = string.join(StoredFileLayout, separator_column)

                with open('%s/StatisticalData/TotalPSData_%s.txt' % (
                    self.storage_options['DATADIR'], StringArgument_other), 'w') as f:
                    for x in zip(*StoredStatisticalData):
                        f.write("%s\n" % (StoredFileLayout).format(*x))

                f.close()

        # Move all the information into subdirectories.
        if self.use_lightcone:

            if self.storage_options['KEEP_GLOBAL_DATA'] or self.FlagOptions.PRINT_FILES:
                LightconePSFilename = 'delTps_lightcone_filenames_%s.txt' % (StringArgument_other)
                filename = open('%s' % (LightconePSFilename), 'r')
                LightconePS = [line.rstrip('\n') for line in filename]

            if self.storage_options['KEEP_ALL_DATA']:
                os.rename(LightconePSFilename,
                          "%s/StatisticalData/%s" % (self.storage_options['DATADIR'], LightconePSFilename))

            if self.FlagOptions.PRINT_FILES:
                os.remove(LightconePSFilename)

            # Removal of the individual light cone files is done here as in principle these can exceed the
            # number of mock observations provided
            if self.storage_options['KEEP_ALL_DATA']:
                for i in range(len(LightconePS)):
                    os.rename(LightconePS[i],
                              "%s/StatisticalData/%s" % (self.storage_options['DATADIR'], LightconePS[i]))

            if self.FlagOptions.PRINT_FILES:
                for i in range(len(LightconePS)):
                    os.remove(LightconePS[i])
        else:

            if self.FlagOptions.PRINT_FILES:
                os.remove("delTps_estimate_%s_*" % StringArgument_other)
                os.remove("NeutralFraction_%s_*" % StringArgument_other)

        if self.FlagOptions.STORE_DATA:
            avefile = "AveData_%s.txt" % StringArgument_other
            if self.storage_options['KEEP_ALL_DATA']:
                os.rename(avefile, "%s/AveData/%s" % (self.storage_options['DATADIR'], avefile))

            if self.FlagOptions.PRINT_FILES:
                os.remove(avefile)

        walkerfile = "Walker_%s.txt" % StringArgument_other
        walkercosmofile = "WalkerCosmology_%s.txt" % StringArgument_other
        if self.storage_options['KEEP_ALL_DATA']:
            os.rename(walkerfile, "%s/WalkerData/%s" % (self.storage_options['DATADIR'], walkerfile))
            os.rename(walkercosmofile, "%s/WalkerData/%s" % (self.storage_options['DATADIR'], walkercosmofile))

        if self.FlagOptions.PRINT_FILES:
            os.remove(walkerfile)
            os.remove(walkercosmofile)
