parameter_scan = [10000, 1000, 100, 10, 1]
time = [50, 50, 100, 100, 100]
total_landscapes = []
fig, gs, axes = generate_axes_pad(nrows=len(parameter_scan), ncols=2,
                                  v_pad=0.4, h_pad=0.2, figsize=(12, 12))
for MC_interval in parameter_scan:
    start = timer()
    print('MC_interval = {}'.format(MC_interval))   
    timesteps = time[parameter_scan.index(MC_interval)]*MC_interval
    walker = empty((total_timesteps/MC_interval, timesteps))
    net_flux = empty((total_timesteps/MC_interval, timesteps))
    steps_executed = []
    while sum(steps_executed) < total_timesteps:
        if i == 0:
            position = pi
        else:
            position = walker[-1][-1]
        this_run = simulate(position, dx, D, kT, dt, shift, forces, energy, i)
        walker[i,0:len(this_run[0])] = this_run[0]
        net_flux[i,0:len(this_run[2])] = this_run[2]
        steps_executed.append(this_run[3])
        i += 1

    simulation_time = timer() - start
    total_landscapes.append(i)
    print('Simulation took {} seconds'.format(simulation_time))
    print('Total landscapes sampled: {}'.format(i+1))

    ################################
    # CONDENSE THE DATA
    ################################

    positions = [walker[i][walker[i]!=0] for i in range(len(walker))]
    # steps = [jumps[i][jumps[i]!=0] for i in range(len(jumps))]
    # Don't have to renormalized net_flux, because extra zeroes are okay.
    fluxes = [sum(net_flux[i]) for i in range(len(net_flux))]
    print('Overall net flux = {}'.format(sum(fluxes)))
    print('Ratio of steps on energy landscape A to B = {}'.
          format(float(len(hstack(positions[0::2]))) /
                 float(len(hstack(positions[1::2])))))
    print('Monte Carlo acceptance ratio = 1 out of every {} attempts'.
          format(mean(steps_executed)/MC_interval))
            

    j = parameter_scan.index(MC_interval)
    c = clrs[j % 9]
    ax = axes[j][0]
    counts, edges = histogram(hstack(positions[0::2]), range=(min(q), max(q)),
                              bins=len(q), density=True)
    mids = (edges[1:] + edges[:-1]) / 2.
    ax.bar(mids, counts, color=c, edgecolor='none',
           width=mids[1] - mids[0], label='{} steps'.format(len(hstack(positions[0::2]))))
    ax.plot(q, pdf[0], color='k', lw=2)
    ax.set_title('Net flux = {}'.format(sum(fluxes)))
    ax.legend()

    ax = axes[j][1]
    # i = 0, then walker[0], walker[2], walker[4]
    # i = 1, then walker[1], walker[3], walker[5]
    counts, edges = histogram(hstack(positions[1::2]), range=(min(q), max(q)),
                              bins=len(q), density=True)
    mids = (edges[1:] + edges[:-1]) / 2.
    ax.bar(mids, counts, color=c, edgecolor='none',
           width=mids[1] - mids[0], label='{} steps'.format(len(hstack(positions[1::2]))))
    ax.plot(q, pdf[1], color='k', lw=2)
    ax.set_title('Net flux = {}, MC interval = {}, landscapes = {}'.format(sum(fluxes), MC_interval, total_landscapes[j]))
    ax.legend()

    ##########################
    # CLEANUP
    ##########################
    del walker
    del net_flux
    del positions

plt.show()
