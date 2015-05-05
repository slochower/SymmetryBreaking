def simulateBDMC(x, dx, D, kT, dt, shift, forces, energies, i, steps_on_this_landscape, **kwargs):

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
