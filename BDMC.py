def initialize_BDMC():

    energies = [energy(q, shift[i]) for i in range(len(shift))]
    if flashing:
        energies[1] = array([mean(energies[0])]*len(energies[0]))
    forces = [[force(k, dx, energies[i]) for k in range(len(q)-1)]
              for i in range(len(energies))]
    boltzmann = [exp(-energies[i]/kT) for i in range(len(energies))]
    pdf = [boltzmann[i] / sum(boltzmann[i]*dx) for i in range(len(boltzmann))]

    walker, net_flux = zeros((2, total_timesteps)), zeros((2, total_timesteps))
    steps_executed, landscapes_sampled, start = 0, 0, 0.0
    steps_on_A, steps_on_B = 0, 0
    return(energies, forces, boltzmann, pdf, walker, net_flux, steps_executed, landscapes_sampled, start, steps_on_A, steps_on_B)


def simulate_BDMC(x, dx, D, kT, dt, shift, forces, energies, i, steps_on_this_landscape, **kwargs):

    positions = []
    fluxes = []
    # Crossings of the flux barrier
    # -1: left crossing of the boundary
    #  0: no crossing this step
    # +1: right crossing of the boundary
    flux_point = pi

    if (i % 2) == 0: 
        state = 0
    else:
        state = 1

    if debug:
        print('\nStarting walker at x = {}'.format(x))
        print('Running for a maximum of {} steps on this landscape'.format(steps_on_this_landscape))
        print('Recording position = {}'.format(x))
    if debug_alternative:
        print('({}, {})'.format(state,x))
    # Record the initial position of the walker at t = 0 and set
    # the flux at this time to be 0, by definition.
    positions.append(x)
    fluxes.append(0)
    # Since we have recorded the position and flux, count that as a 
    # time step.
    t = 1
    # Each iteration through the loop adds two timesteps: one for MC and
    # one for BD, so if there is only 1 step remaining, then quit early.
    quit_early = False
    if steps_on_this_landscape == 1:
        quit_early = True
    while t < steps_on_this_landscape-1:
        #######################################
        # BROWNIAN WALK
        #######################################
        g = random.normal(loc=0.0, scale=sqrt(2 * D * dt))
        F = force_lookup(q[:-1], forces[state], x)
        new_x = x + (D / kT) * F * dt + g

        if (abs(new_x - x) > max(q)):
            raise Exception('Jumps are too big.')

        #######################################
        # MEASURE FLUX
        #######################################
        if (x < flux_point) and (new_x >= flux_point):
            swap = 1
        elif (x > flux_point) and (new_x <= flux_point):
            swap = -1
        else:
            swap = 0

        ########################################
        # KEEP TRACK OF PBC CROSSINGS
        ########################################
        if new_x > max(q):
            new_x = min(q) + (new_x - max(q))
        elif new_x < min(q):
            new_x = max(q) - abs(new_x - min(q))
        if debug:
            print('t = {}, landscape = {}, x = {}, new_x = {}, status = {}'.
                  format(t, state, x, new_x, str('BD')))
            print('Recording position = {}'.format(new_x))
        if debug_alternative:
            print('BD')
            print('({}, {})'.format(state,new_x))
        ####################
        # RECORD KEEPING
        ####################
        t += 1
        fluxes.append(swap)
        positions.append(new_x)
        x = new_x
        if quit_early: break
        #######################################
        # MONTE CARLO CHECK TO STEP ORTHOGONAL
        #######################################
        if (MC is False):
            raise Execption('Brownian dynamics without MC is not implemented.')

        if (MC is True and ((t+1) % MC_interval) == 0):
            E_transition = energy_lookup(q, energies[(state+1) % 2], x)
            E_now = energy_lookup(q, energies[state], x)
            delta = E_transition - E_now
            # If delta < 0, step to next landscape
            if (delta < 0):
                if debug:
                    print('t = {}, landscape = {}, x = {}, status = {}'.
                          format(t, state, x, str('MC ACCEPT (delta E)')))
                if debug_alternative:
                    print('MC accept')
                #    print('Recording position = {}'.format(x))
                #positions.append(x)
                #fluxes.append(0)
                #t += 1
                break
            # If delta !< 0, compute p_accept and pick a random number.
            p_accept = exp(-delta/kT)
            r = random.random()
            if (p_accept > r):
                if debug:
                    print('t = {}, landscape = {}, x = {}, status = {}'.
                          format(t, state, x, str('MC ACCEPT (p_accept > r)')))
                if debug_alternative:
                    print('MC accept')

                #    print('Recording position = {}'.format(x))
                #positions.append(x)
                #fluxes.append(0)
                #t += 1
                break
            # If p_accept !> r, append the current position and 
            # then take a Brownian dynamics step.
            if debug:
                print('t = {}, landscape = {}, x = {}, status = {}'.
                      format(t, state, x, str('MC FAIL')))
                print('Recording position = {}'.format(x))
            if debug_alternative:
                print('MC fail')
                print('({}, {})'.format(state,x))

            t += 1
            positions.append(x)
            fluxes.append(0)
            pass

    return(positions, fluxes)



def simulate_BDMC_both_surfaces(x, dx, D, kT, dt, shift, forces, energies, i, steps_on_this_landscape, **kwargs):

    positions = []
    fluxes = []
    # Crossings of the flux barrier
    # -1: left crossing of the boundary
    #  0: no crossing this step
    # +1: right crossing of the boundary
    MC_side_attempt = 0
    MC_orthogonal_attempt = 0

    flux_point = pi

    if (i % 2) == 0: 
        state = 0
    else:
        state = 1

    if debug:
        print('\nStarting walker at x = {}'.format(x))
        print('Running for a maximum of {} steps on this landscape'.format(steps_on_this_landscape))
        print('Recording position = {}'.format(x))
    if debug_alternative:
        print('({}, {})'.format(state,x))
    # Record the initial position of the walker at t = 0 and set
    # the flux at this time to be 0, by definition.
    positions.append(x)
    fluxes.append(0)
    # Since we have recorded the position and flux, count that as a 
    # time step.
    t = 1
    # Each iteration through the loop adds two timesteps: one for MC and
    # one for BD, so if there is only 1 step remaining, then quit early.
    quit_early = False
    if steps_on_this_landscape == 1:
        quit_early = True
    while t < steps_on_this_landscape-1:
        #######################################
        # BROWNIAN WALK
        #######################################
        g = random.normal(loc=0.0, scale=sqrt(2 * D * dt))
        F = force_lookup(q[:-1], forces[state], x)
        new_x = x + (D / kT) * F * dt + g

        if (abs(new_x - x) > max(q)):
            raise Exception('Jumps are too big.')

        #######################################
        # MEASURE FLUX
        #######################################
        if (x < flux_point) and (new_x >= flux_point):
            swap = 1
        elif (x > flux_point) and (new_x <= flux_point):
            swap = -1
        else:
            swap = 0

        ########################################
        # KEEP TRACK OF PBC CROSSINGS
        ########################################
        if new_x > max(q):
            new_x = min(q) + (new_x - max(q))
        elif new_x < min(q):
            new_x = max(q) - abs(new_x - min(q))
        if debug:
            print('t = {}, landscape = {}, x = {}, new_x = {}, status = {}'.
                  format(t, state, x, new_x, str('BD')))
            print('Recording position = {}'.format(new_x))
        if debug_alternative:
            print('BD')
            print('({}, {})'.format(state,new_x))
        ####################
        # RECORD KEEPING
        ####################
        t += 1
        fluxes.append(swap)
        positions.append(new_x)
        x = new_x
        if quit_early: break
        #######################################
        # MONTE CARLO CHECK TO STEP ORTHOGONAL
        #######################################
        if (MC is False):
            raise Execption('Brownian dynamics without MC is not implemented.')

        if (MC is True and ((t+1) % MC_interval) == 0):
            ###########################################################
            # WALK ORTHOGONAL?
            ###########################################################
            r = random.random()
            if r > 0.5:
                E_transition = energy_lookup(q, energies[(state+1) % 2], x)
                E_now = energy_lookup(q, energies[state], x)
                delta = E_transition - E_now
                # If delta < 0, step to next landscape
                if (delta < 0):
                    if debug:
                        print('t = {}, landscape = {}, x = {}, status = {}'.
                              format(t, state, x, str('MC ACCEPT (delta E)')))
                    if debug_alternative:
                        print('MC accept')
                    #    print('Recording position = {}'.format(x))
                    #positions.append(x)
                    #fluxes.append(0)
                    #t += 1
                    break
                # If delta !< 0, compute p_accept and pick a random number.
                p_accept = exp(-delta/kT)
                r = random.random()
                if (p_accept > r):
                    if debug:
                        print('t = {}, landscape = {}, x = {}, status = {}'.
                              format(t, state, x, str('MC ACCEPT (p_accept > r)')))
                    if debug_alternative:
                        print('MC accept')

                    #    print('Recording position = {}'.format(x))
                    #positions.append(x)
                    #fluxes.append(0)
                    #t += 1
                    break
                # If p_accept !> r, append the current position and 
                # then take a Brownian dynamics step.
                if debug:
                    print('t = {}, landscape = {}, x = {}, status = {}'.
                          format(t, state, x, str('MC FAIL')))
                    print('Recording position = {}'.format(x))
                if debug_alternative:
                    print('MC fail')
                    print('({}, {})'.format(state,x))

                t += 1
                positions.append(x)
                fluxes.append(0)
                MC_orthogonal_attempt += 1
                pass
            else:
                #####################################################
                # TRY MC ON CURRENT LANDSCAPE
                #####################################################

                b =  0.2
                a = -0.2
                new_x = ((b - a) * random.random() + a) + x
                if new_x > max(q):
                    new_x = min(q) + (new_x - max(q))
                elif new_x < min(q):
                    new_x = max(q) - abs(new_x - min(q))

                if debug_alternative:
                    print('Trial step = {}'.format(new_x))
                E_step = energy_lookup(q, energies[state], new_x)
                E_now = energy_lookup(q, energies[state], x)
                delta = E_step - E_now
                if (delta < 0):
                    positions.append(new_x)
                    if debug_alternative:
                        print('MC forward step accepted (delta E)')

                else:
                    p_accept = exp(-delta/kT)
                    r = random.random()
                    if (p_accept > r):
                        positions.append(new_x)
                        if debug_alternative:
                            print('MC forward step accepted (random)')
                    else:
                        new_x = x
                        positions.append(new_x)

                if (abs(new_x - x) > max(q)):
                    raise Exception('Jumps are too big.')

                #######################################
                # MEASURE FLUX
                #######################################
                if (x < flux_point) and (new_x >= flux_point):
                    swap = 1
                elif (x > flux_point) and (new_x <= flux_point):
                    swap = -1
                else:
                    swap = 0

                if debug:
                    print('t = {}, landscape = {}, x = {}, new_x = {}, status = {}'.
                          format(t, state, x, new_x, str('BD')))
                    print('Recording position = {}'.format(new_x))
                if debug_alternative:
                    print('({}, {})'.format(state,new_x))
                ####################
                # RECORD KEEPING
                ####################
                t += 1
                fluxes.append(swap)
                x = new_x
                MC_side_attempt += 1


    return(positions, fluxes, MC_orthogonal_attempt, MC_side_attempt)
