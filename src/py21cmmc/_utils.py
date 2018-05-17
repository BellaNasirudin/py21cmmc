
def write_walker_file(flag_options, astro_params, random_ids= ('3.000000', '3.000000')):
    separator = " "
    separator_other = "_"
    seq = []

    #TODO: this is probably better as a hash or something. Add the random thread ID
    seq.append("%s" % (random_ids[0]))
    # Add the second ID
    seq.append("%s" % (random_ids[1]))


    StringArgument_other = separator_other.join(seq)

    # Add number of redshifts
    # If using the light-cone version of the code, don't need to set a redshift
    seq.append("%d" % (flag_options.N_USER_REDSHIFTS))
    # Add light cone flag
    seq.append("%d" % (flag_options.USE_LIGHTCONE))
    # If power-law dependence on ionising efficiency is allowed. Add the flag here (support not yet included)
    seq.append("0")
    # Add redshift for Ts.c calculation
    seq.append("%g" % (flag_options.REDSHIFT))

    StringArgument = separator.join(seq)

    with open("Walker_%s.txt" % (StringArgument_other), "w") as f:
        f.write("FLAGS")
        for key in flag_options._fields_:
            if key[0] not in ['N_USER_REDSHIFTS', 'USE_LIGHTCONE', 'REDSHIFT']:
                f.write("    %s"%(getattr(flag_options, key[0])))
        f.write("\n")

        for key in astro_params._fields_:
            f.write("%s    %f\n"%(key[0],getattr(astro_params, key[0])))

        if hasattr(flag_options.redshifts, "__len__"):
            for z in flag_options.redshifts:
                f.write("CO-EVAL-Z    %f\n" % (float(z)))


def write_walker_cosmology_file( flag_options, cosmo_params, random_ids):
    separator_other = "_"
    seq = []

    #TODO: this is probably better as a hash or something. Add the random thread ID
    seq.append("%s" % (random_ids[0]))
    # Add the second ID
    seq.append("%s" % (random_ids[1]))


    StringArgument_other = separator_other.join(seq)

    # Add number of redshifts
    # If using the light-cone version of the code, don't need to set a redshift
    seq.append("%d" % (flag_options.N_USER_REDSHIFTS))
    # Add light cone flag
    seq.append("%d" % (flag_options.USE_LIGHTCONE))
    # If power-law dependence on ionising efficiency is allowed. Add the flag here (support not yet included)
    seq.append("0")
    # Add redshift for Ts.c calculation
    seq.append("%g" % (flag_options.REDSHIFT))

    with open("WalkerCosmology_%s.txt"%StringArgument_other, 'w') as f:
        for k in cosmo_params._fields_:
            f.write("%s    %f\n"%(k[0], getattr(cosmo_params, k[0])))
