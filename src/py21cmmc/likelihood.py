"""
A module containing (base) classes for computing 21cmFAST likelihoods under the context of CosmoHammer.
"""
import numpy as np
from decimal import *
from scipy import interpolate
from scipy.interpolate import InterpolatedUnivariateSpline
from . import wrapper as p21c

import pickle
from os import path

np.seterr(invalid='ignore', divide='ignore')

try:
    from powerbox.tools import get_power
    from powerbox import dft
    from astropy.cosmology import Planck15

    HAVE_PB_AP = True
except ImportError:
    HAVE_PB_AP = False

TWOPLACES = Decimal(10) ** -2  # same as Decimal('0.01')
FOURPLACES = Decimal(10) ** -4  # same as Decimal('0.0001')
SIXPLACES = Decimal(10) ** -6  # same as Decimal('0.000001')


class LikelihoodBase:
    def __init__(self, box_dim={}, flag_options={}, astro_params={}, cosmo_params={}):
        # Save the following only as dictionaries (not full Structure objects), so that we can avoid
        # weird pickling errors.
        self._box_dim = box_dim
        self._flag_options = flag_options
        self._astro_params = astro_params
        self._cosmo_params = cosmo_params

        self.use_lightcone = self.flag_options().USE_LIGHTCONE
        self.redshift = self.flag_options().redshifts

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
        return p21c.FlagOptionStruct(self.astro_params().Z_HEAT_MAX, self.astro_params().ZPRIME_STEP_FACTOR,
                                     **dct)

    def computeLikelihood(self, ctx):
        raise NotImplementedError("The Base likelihood should never be used directly!")

    def setup(self):
        pass


class LikelihoodPlanck(LikelihoodBase):
    # Mean and one sigma errors for the Planck constraints
    # The Planck prior is modelled as a Gaussian: tau = 0.058 \pm 0.012 (https://arxiv.org/abs/1605.03507)
    PlanckTau_Mean = 0.058
    PlanckTau_OneSigma = 0.012

    # Simple linear extrapolation of the redshift range provided by the user, to be able to estimate the optical depth
    nZinterp = 15

    # The minimum of the extrapolation is chosen to 5.9, to correspond to the McGreer et al. prior on the IGM neutral fraction.
    # The maximum is chosed to be z = 18., which is arbitrary.
    ZExtrap_min = 5.9
    ZExtrap_max = 20.0

    def computeLikelihood(self, ctx):
        """
        Contribution to the likelihood arising from Planck (2016) (https://arxiv.org/abs/1605.03507)
        """
        # READ_FROM_FILE = ctx.get('flag_options').READ_FROM_FILE
        # PRINT_FILES = ctx.get('FlagOptions').PRINT_FILES

        # Extract relevant info from the context.
        output = ctx.get("output")

        if len(output.redshifts) < 3:
            print(output.redshifts)
            raise ValueError("You cannot use the Planck prior likelihood with less than 3 redshifts")

        # The linear interpolation/extrapolation function, taking as input the redshifts supplied by the user and
        # the corresponding neutral fractions recovered for the specific EoR parameter set
        LinearInterpolationFunction = InterpolatedUnivariateSpline(output.redshifts, output.average_nf, k=1)

        ZExtrapVals = np.zeros(self.nZinterp)
        XHI_ExtrapVals = np.zeros(self.nZinterp)

        for i in range(self.nZinterp):
            ZExtrapVals[i] = self.ZExtrap_min + (self.ZExtrap_max - self.ZExtrap_min) * float(i) / (self.nZinterp - 1)

            XHI_ExtrapVals[i] = LinearInterpolationFunction(ZExtrapVals[i])

            # Ensure that the neutral fraction does not exceed unity, or go negative
            if XHI_ExtrapVals[i] > 1.0:
                XHI_ExtrapVals[i] = 1.0
            if XHI_ExtrapVals[i] < 0.0:
                XHI_ExtrapVals[i] = 0.0

        # Set up the arguments for calculating the estimate of the optical depth. Once again, performed using command line code.
        tau_value = p21c.compute_tau(ZExtrapVals, XHI_ExtrapVals, ctx.get('cosmo_params'))

        # remove the temporary files (this depends on tau being run, so don't move it to _store_data())
        # if self.FlagOptions.PRINT_FILES:
        #     taufile = "Tau_e_%s_%s.txt" % random_ids
        #     if self.storage_options['KEEP_ALL_DATA']:
        #         os.rename(taufile, "%s/TauData/%s" % (self.storage_options['DATADIR'], taufile))
        #     else:
        #         os.remove(taufile)

        # As the likelihood is computed in log space, the addition of the prior is added linearly to the existing chi^2 likelihood
        lnprob = np.square((self.PlanckTau_Mean - tau_value) / (self.PlanckTau_OneSigma))

        return lnprob

        # TODO: not sure what to do about this:
        # it is len(self.AllRedshifts) as the indexing begins at zero


#        nf_vals[len(self.AllRedshifts) + 2] = tau_value


class LikelihoodMcGreer(LikelihoodBase):
    # Mean and one sigma errors for the McGreer et al. constraints
    # Modelled as a flat, unity prior at x_HI <= 0.06, and a one sided Gaussian at x_HI > 0.06
    # ( Gaussian of mean 0.06 and one sigma of 0.05 )
    McGreer_Mean = 0.06
    McGreer_OneSigma = 0.05
    McGreer_Redshift = 5.9

    def computeLikelihood(self, ctx):
        """
        Limit on the IGM neutral fraction at z = 5.9, from dark pixels by I. McGreer et al.
        (2015) (http://adsabs.harvard.edu/abs/2015MNRAS.447..499M)
        """
        lightcone = ctx.get("output")

        if self.McGreer_Redshift in lightcone.redshifts:
            for i in range(len(lightcone.redshifts)):
                if lightcone.redshifts[i] == self.McGreer_Redshift:
                    McGreer_NF = lightcone.average_nf[i]
        elif len(lightcone.redshifts) > 2:
            # The linear interpolation/extrapolation function, taking as input the redshifts supplied by the user and
            # the corresponding neutral fractions recovered for the specific EoR parameter set
            LinearInterpolationFunction = InterpolatedUnivariateSpline(lightcone.redshifts, lightcone.average_nf, k=1)
            McGreer_NF = LinearInterpolationFunction(self.McGreer_Redshift)
        else:
            raise ValueError(
                "You cannot use the McGreer prior likelihood with either less than 3 redshifts or the redshift being directly evaluated.")

        McGreer_NF = np.clip(McGreer_NF, 0, 1)

        lnprob = 0
        if McGreer_NF > 0.06:
            lnprob = np.square((self.McGreer_Mean - McGreer_NF) / (self.McGreer_OneSigma))

        return lnprob


class LikelihoodGreig(LikelihoodBase):
    QSO_Redshift = 7.0842  # The redshift of the QSO

    def __init__(self, *args, **kwargs):

        super().__init__(*args, **kwargs)

    def setup(self):
        with open(path.expanduser(path.join("~", '.py21cmmc', 'PriorData', "NeutralFractionsForPDF.out")),
                  'rb') as handle:
            self.NFValsQSO = pickle.loads(handle.read())

        with open(path.expanduser(path.join("~", '.py21cmmc', 'PriorData', "NeutralFractionPDF_SmallHII.out")),
                  'rb') as handle:
            self.PDFValsQSO = pickle.loads(handle.read())

        # Normalising the PDF to have a peak probability of unity (consistent with how other priors are treated)
        # Ultimately, this step does not matter
        normalisation = np.amax(self.PDFValsQSO)
        self.PDFValsQSO /= normalisation

    def computeLikelihood(self, ctx):
        """
        Constraints on the IGM neutral fraction at z = 7.1 from the IGM damping wing of ULASJ1120+0641
        Greig et al (2016) (http://arxiv.org/abs/1606.00441)
        """

        lightcone = ctx.get("output")

        Redshifts = lightcone.redshifts
        AveNF = lightcone.average_nf

        # Interpolate the QSO damping wing PDF
        spline_QSODampingPDF = interpolate.splrep(self.NFValsQSO, self.PDFValsQSO, s=0)

        if self.QSO_Redshift in Redshifts:

            for i in range(len(Redshifts)):
                if Redshifts[i] == self.QSO_Redshift:
                    NF_QSO = AveNF[i]

        elif len(lightcone.redshifts) > 2:

            # Check the redshift range input by the user to determine whether to interpolate or extrapolate the IGM
            # neutral fraction to the QSO redshift
            if self.QSO_Redshift < np.amin(Redshifts):
                # The QSO redshift is outside the range set by the user. Need to extrapolate the reionisation history
                # to obtain the neutral fraction at the QSO redshift

                # The linear interpolation/extrapolation function, taking as input the redshifts supplied by the user
                # and the corresponding neutral fractions recovered for the specific EoR parameter set
                LinearInterpolationFunction = InterpolatedUnivariateSpline(Redshifts, AveNF, k=1)

                NF_QSO = LinearInterpolationFunction(self.QSO_Redshift)

            else:
                # The QSO redshift is within the range set by the user. Can interpolate the reionisation history to
                # obtain the neutral fraction at the QSO redshift
                if lightcone.params.n_redshifts == 3:
                    spline_reionisationhistory = interpolate.splrep(Redshifts, AveNF, k=2, s=0)
                else:
                    spline_reionisationhistory = interpolate.splrep(Redshifts, AveNF, s=0)

                NF_QSO = interpolate.splev(self.QSO_Redshift, spline_reionisationhistory, der=0)

        else:
            raise ValueError(
                "You cannot use the Greig prior likelihood with either less than 3 redshifts or the redshift being directly evaluated.")

        # Ensure that the neutral fraction does not exceed unity, or go negative
        NF_QSO = np.clip(NF_QSO, 0, 1)

        QSO_Prob = interpolate.splev(NF_QSO, spline_QSODampingPDF, der=0)

        # Interpolating the PDF from the QSO damping wing might cause small negative values at the edges (i.e. x_HI ~ 0 or ~1)
        # In case it is zero, or negative, set it to a very small non zero number (we take the log of this value, it cannot be zero)
        if QSO_Prob <= 0.0:
            QSO_Prob = 0.000006

        # We work with the log-likelihood, therefore convert the IGM Damping wing PDF to log space
        QSO_Prob = -2. * np.log(QSO_Prob)

        lnprob = QSO_Prob
        return lnprob


class LikelihoodGlobal(LikelihoodBase):

    def __init__(self, FIXED_ERROR=False,
                 model_name="FaintGalaxies", mock_dir=None,
                 fixed_global_error=10.0, fixed_global_bandwidth=4.0, FrequencyMin=40.,
                 FrequencyMax=200, *args, **kwargs):

        # Run the LikelihoodBase init.
        super().__init__(*args, **kwargs)

        self.FIXED_ERROR = FIXED_ERROR

        self.model_name = model_name
        self.mock_dir = mock_dir or path.expanduser(path.join("~", '.py21cmmc'))

        self.fixed_global_error = fixed_global_error
        self.fixed_global_bandwidth = fixed_global_bandwidth
        self.FrequencyMin = FrequencyMin
        self.FrequencyMax = FrequencyMax

        self.obs_filename = path.join(self.mock_dir, "MockData", self.model_name, "GlobalSignal",
                                      self.model_name + "_GlobalSignal.txt")
        self.obs_error_filename = path.join(self.mock_dir, 'NoiseData', self.model_name, "GlobalSignal",
                                            'TotalError_%s_GlobalSignal_ConstantError_1000hr.txt' % self.model_name)

    def setup(self):
        """
        Contains any setup specific to this likelihood, that should only be performed once. Must save variables
        to the class.
        """

        # Read in the mock 21cm PS observation. Read in both k and the dimensionless PS.
        self.k_values = []
        self.PS_values = []

        mock = np.loadtxt(self.obs_filename, usecols=(0, 2))
        self.k_values.append(mock[:, 0])
        self.PS_values.append(mock[:, 1])

        self.Error_k_values = []
        self.PS_Error = []

        if not self.FIXED_ERROR:
            errs = np.loadtxt(self.obs_error_filename, usecols=(0, 1))

            self.Error_k_values.append(errs[:, 0])
            self.PS_Error.append(errs[:, 1])

        self.Error_k_values = np.array(self.Error_k_values)
        self.PS_Error = np.array(self.PS_Error)

    def computeLikelihood(self, ctx):
        """
        Compute the likelihood, given the lightcone output from 21cmFAST.
        """
        lightcone = ctx.get("output")

        # Get some useful variables out of the Lightcone box
        NumRedshifts = len(lightcone.redshifts)
        Redshifts = lightcone.redshifts
        AveTb = lightcone.average_Tb

        total_sum = 0

        # Converting the redshifts to frequencies for the interpolation (must be in increasing order, it is by default redshift which is decreasing)
        FrequencyValues_mock = np.zeros(len(self.k_values[0]))
        FrequencyValues_model = np.zeros(NumRedshifts)

        # Shouldn't need two, as they should be the same sampling. However, just done it for now
        for j in range(len(self.k_values[0])):
            FrequencyValues_mock[j] = ((2.99792e8) / (.2112 * (1. + self.k_values[0][j]))) / (1e6)

        for j in range(NumRedshifts):
            FrequencyValues_model[j] = ((2.99792e8) / (.2112 * (1. + Redshifts[j]))) / (1e6)

        splined_mock = interpolate.splrep(FrequencyValues_mock, self.PS_values[0], s=0)
        splined_model = interpolate.splrep(FrequencyValues_model, AveTb, s=0)

        FrequencyMin = self.FrequencyMin
        FrequencyMax = self.FrequencyMax

        if self.FIXED_ERROR:
            ErrorOnGlobal = self.fixed_global_error
            Bandwidth = self.fixed_global_bandwidth

            FrequencyBins = int(np.floor((FrequencyMax - FrequencyMin) / Bandwidth)) + 1

            for j in range(FrequencyBins):
                FrequencyVal = FrequencyMin + Bandwidth * j

                MockPS_val = interpolate.splev(FrequencyVal, splined_mock, der=0)

                ModelPS_val = interpolate.splev(FrequencyVal, splined_model, der=0)

                total_sum += np.square((MockPS_val - ModelPS_val) / ErrorOnGlobal)

        else:

            for j in range(len(self.Error_k_values[0])):

                FrequencyVal = ((2.99792e8) / (.2112 * (1. + self.Error_k_values[0][j]))) / (1e6)

                if FrequencyVal >= FrequencyMin and FrequencyVal <= FrequencyMax:
                    MockPS_val = interpolate.splev(FrequencyVal, splined_mock, der=0)

                    ModelPS_val = interpolate.splev(FrequencyVal, splined_model, der=0)

                    total_sum += np.square((MockPS_val - ModelPS_val) / self.PS_Error[0][j])

        return -0.5 * total_sum  # , nf_vals


class Likelihood1DPowerMultiZ(LikelihoodBase):

    def __init__(self, data_redshifts=None, NSplinePoints=8, Foreground_cut=0.15, Shot_Noise_cut=1.0,
                 ModUncert=0.2, log_sampling=False, model_name="FaintGalaxies", mock_dir=None,
                 telescope="HERA331", duration='1000hr',
                 *args, **kwargs):

        # Run the LikelihoodBase init.
        super().__init__(*args, **kwargs)

        self.Foreground_cut = Foreground_cut
        self.Shot_Noise_cut = Shot_Noise_cut
        self.NSplinePoints = NSplinePoints
        self.ModUncert = ModUncert
        self.data_redshifts = data_redshifts
        self.log_sampling = log_sampling

        self.telescope = telescope
        self.duration = duration
        self.model_name = model_name
        self.mock_dir = mock_dir or path.expanduser(path.join("~", '.py21cmmc'))

    def setup(self):
        """
        Contains any setup specific to this likelihood, that should only be performed once. Must save variables
        to the class.
        """
        if self.log_sampling:
            kSplineMin = np.log10(self.Foreground_cut)
            kSplineMax = np.log10(self.Shot_Noise_cut)
        else:
            kSplineMin = self.Foreground_cut
            kSplineMax = self.Shot_Noise_cut

        kSpline = np.zeros(self.NSplinePoints)

        # TODO: probably should be np.linspace.
        for j in range(self.NSplinePoints):
            kSpline[j] = kSplineMin + (kSplineMax - kSplineMin) * float(j) / (self.NSplinePoints - 1)

        if self.log_sampling:
            self.kSpline = 10 ** kSpline
        else:
            self.kSpline = kSpline

        self.k_values, self.PS_values, self.Error_k_values, self.PS_Error = self.define_data()

    def define_data(self):
        """
        An over-rideable method which should return the k, P, and error on power spectrum at each redshift. Nominally,
        reads these in from files.

        Returns
        -------
        k_values, PS_values, Error_k_values, PS_Error
            Must return these values.
        """

        if self.use_lightcone:
            self.obs_filename = path.join(self.mock_dir, 'MockData',
                                          'LightCone21cmPS_%s_600Mpc_400.txt' % self.model_name)
            self.obs_error_filename = path.join(self.mock_dir, 'NoiseData',
                                                'LightCone21cmPS_Error_%s_%s_%s_600Mpc_400.txt' % (
                                                self.model_name, self.telescope, self.duration))

        else:
            if not self.data_redshifts:
                raise ValueError("If not using a lightcone, you must pass at least some data redshifts")

            if any([z not in self.redshift for z in self.data_redshifts]):
                raise ValueError("One or more data redshifts were not in the computed redshifts %s %s." % (
                self.data_redshifts, self.redshift))

            self.obs_filename = path.join(self.mock_dir, 'MockData', self.model_name, "Co-Eval",
                                          'MockObs_%s_PS_200Mpc_' % self.model_name)
            self.obs_error_filename = path.join(self.mock_dir, 'NoiseData', self.model_name, "Co-Eval",
                                                'TotalError_%s_PS_200Mpc.txt' % self.telescope)

        if not path.exists(self.obs_filename) or not path.exists(self.obs_error_filename):
            raise ValueError("Those mock observations and/or noise files do not exist: %s %s" % (
            self.obs_filename, self.obs_error_filename))

        # Read in the mock 21cm PS observation. Read in both k and the dimensionless PS.
        # These are needed for performing the chi^2 statistic for the likelihood. NOTE: To calculate the likelihood
        # statistic a spline is performed for each of the mock PS, simulated PS and the Error PS
        k_values = []
        PS_values = []

        if self.use_lightcone:
            # Note here, we are populating the list 'Redshift' with the filenames. The length of this is needed for
            # ensuring the correct number of 21cm PS are used for the likelihood. Re-using the same list filename
            # means less conditions further down this script. The likelihood correctly accounts for it with the
            # 'use_lightcone' flag.
            with open(self.obs_filename, 'r') as f:
                subfiles = [line.rstrip('\n') for line in f]

            for fl in subfiles:
                mock = np.loadtxt('%s/%s' % (path.dirname(self.obs_filename), fl), usecols=(0, 1))
                k_values.append(mock[:, 0])
                PS_values.append(mock[:, 1])

        else:

            ### NOTE ###
            # If Include_Ts_fluc is set, the user must ensure that the co-eval redshift to be sampled (set by the
            # Redshift list above) is to be sampled by the code.

            for i, z in enumerate(self.data_redshifts):
                mock = np.loadtxt(self.obs_filename + "%s" % z, usecols=(0, 1))

                k_values.append(mock[:, 0])
                PS_values.append(mock[:, 1])

        k_values = np.array(k_values)
        PS_values = np.array(PS_values)

        ###### Read in the data for the telescope sensitivites ######
        Error_k_values = []
        PS_Error = []

        # Total noise sensitivity as computed from 21cmSense.
        if self.use_lightcone:
            with open(self.obs_error_filename, 'r') as f:
                LightConeErrors = [line.rstrip('\n') for line in f]

            # Use LightConeSnapShots here to ensure it crashes if the number of error files is less than the number or observations
            for i in range(len(subfiles)):
                errs = np.loadtxt('%s/%s' % (path.dirname(self.obs_error_filename), LightConeErrors[i]), usecols=(0, 1))

                Error_k_values.append(errs[:, 0])
                PS_Error.append(errs[:, 1])

        else:

            for i in range(len(self.data_redshifts)):
                errs = np.loadtxt(self.obs_error_filename, usecols=(0, 1))
                Error_k_values.append(errs[:, 0])
                PS_Error.append(errs[:, 1])

        Error_k_values = np.array(Error_k_values)
        PS_Error = np.array(PS_Error)
        return k_values, PS_values, Error_k_values, PS_Error

    def computeLikelihood(self, ctx):
        """
        Compute the likelihood, given the lightcone output from 21cmFAST.
        """
        lightcone = ctx.get("output")

        # Get some useful variables out of the Lightcone box
        PS_Data = lightcone.power_spectrum
        k_Data = lightcone.k

        total_sum = 0

        print(lightcone.power_spectrum, lightcone.k)
        # Note here that the usage of len(redshift) uses the number of mock lightcone 21cm PS if use_lightcone was set to True.
        for i, z in enumerate(self.data_redshifts):

            if not self.use_lightcone:
                redshift_index = np.where(lightcone.redshifts == z)[0][0]
            else:
                redshift_index = i

            splined_mock = interpolate.splrep(self.k_values[i], np.log10(self.PS_values[i]), s=0)
            splined_error = interpolate.splrep(self.Error_k_values[i], np.log10(self.PS_Error[i]), s=0)

            splined_model = interpolate.splrep(k_Data, np.log10(PS_Data[redshift_index]), s=0)

            # Interpolating the mock and error PS in log space
            for j in range(self.NSplinePoints):

                MockPS_val = 10 ** (interpolate.splev(self.kSpline[j], splined_mock, der=0))
                ErrorPS_val = 10 ** (interpolate.splev(self.kSpline[j], splined_error, der=0))

                ModelPS_val = 10 ** (interpolate.splev(self.kSpline[j], splined_model, der=0))

                # Check if there are any nan values for the 21cm PS
                # A nan value implies a IGM neutral fraction of zero, that is, reionisation has completed and thus no 21cm signal
                # Set the value of the 21cm PS to zero. Which results in the largest available difference (i.e. if you expect a signal
                # (i.e. non zero mock 21cm PS) but have no signal from the sampled model, then want a large difference for the
                # chi-squared likelihood).
                if np.isnan(ModelPS_val) == True:
                    ModelPS_val = 0.0

                if np.isnan(MockPS_val) == True:
                    MockPS_val = 0.0

                total_sum += np.square((MockPS_val - ModelPS_val) / (
                    np.sqrt(ErrorPS_val ** 2. + (self.ModUncert * ModelPS_val) ** 2.)))

        return -0.5 * total_sum  # , nf_vals


class Likelihood1DPowerNoErrors(Likelihood1DPowerMultiZ):
    """
    A simple likelihood model that generates "data" as a simple power spectrum from fiducial parameters,
    and applies no noise. Use for testing.
    """

    def define_data(self):
        output = p21c.run_21cmfast(self._flag_options['redshifts'], self._box_dim, self._flag_options,
                                   self._astro_params, self._cosmo_params)[0]

        print(output.power_spectrum, output.k)
        nz = len(self._flag_options['redshifts'])
        return np.repeat(output.k, nz).reshape((len(output.k), nz)).T, output.power_spectrum, np.repeat(output.k,
                                                                                                        nz).reshape(
            (len(output.k), nz)).T, np.ones_like(output.power_spectrum)


class Likelihood1DPowerLightconeNoErrors(LikelihoodBase):
    def __init__(self, datafile, n_psbins=None, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.datafile = datafile
        self.n_psbins = n_psbins

    def setup(self):
        if not HAVE_PB_AP:
            raise NotImplementedError("You need to install powerbox and astropy to use this class")

        if not self.flag_options().USE_LIGHTCONE:
            raise ValueError("You need to use a lightcone for this Likelihood module")

        data = np.genfromtxt(self.datafile)
        self.k_data = data[:, 0]
        self.p_data = data[:, 1]

    def computeLikelihood(self, ctx):
        output = ctx.get("output")
        p, k, var = self.compute_power(output, n_psbins=self.n_psbins)

        return - 0.5 * np.sum((p - self.p_data) ** 2/ var)

    @staticmethod
    def compute_power(lightcone, n_psbins=None):
        p, k, var = get_power(
            lightcone.lightcone_box,
            [lightcone.box_len, lightcone.box_len, lightcone.box_len_lightcone],
            bins=n_psbins, bin_ave=False, get_variance=True
        )

        k = (k[1:] + k[:-1])/2
        return p, k, var

    def simulate_data(self, save_lightcone=False):
        output = p21c.run_21cmfast(self._flag_options['redshifts'], self._box_dim, self._flag_options,
                                   self._astro_params, self._cosmo_params)[0]

        if save_lightcone:
            with open(path.join(path.dirname(self.datafile), "lightcone_data.pkl"), 'wb') as f:
                pickle.dump(output, f)

        p, k, var = self.compute_power(output, n_psbins=self.n_psbins)

        np.savetxt(self.datafile, np.array([k, p]).T)
