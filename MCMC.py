def initialize_MCMC():

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


def simulate_MCMC(x, dx, D, kT, dt, shift, forces, energies, i, steps_on_this_landscape, **kwargs):
    positions = []
    fluxes = []
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
    positions.append(x)
    fluxes.append(0)
    t = 1
    if steps_on_this_landscape == 1:
        quit_early = True
    else:
        quit_early = False

    while t < steps_on_this_landscape-1:
        ##################################
        # MONTE CARLO DYNAMICS
        ##################################
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

        if quit_early: break
        #######################################
        # MONTE CARLO CHECK TO STEP ORTHOGONAL
        #######################################
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
                    print('MC side step accepted (delta E)')
                break
            # If delta !< 0, compute p_accept and pick a random number.
            p_accept = exp(-delta/kT)
            r = random.random()
            if (p_accept > r):
                if debug:
                    print('t = {}, landscape = {}, x = {}, status = {}'.
                          format(t, state, x, str('MC ACCEPT (p_accept > r)')))
                if debug_alternative:
                    print('MC side step accepted (random)')
                break
            # If p_accept !> r, append the current position and 
            # then take a Brownian dynamics step.
            if debug:
                print('t = {}, landscape = {}, x = {}, status = {}'.
                      format(t, state, x, str('MC FAIL')))
                print('Recording position = {}'.format(x))
            if debug_alternative:
                print('MC side step rejected')

            t += 1
            positions.append(x)
            fluxes.append(0)
            pass
    return(positions, fluxes)