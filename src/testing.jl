p = Particle()

bt = billiard_sinai(setting="periodic")

ct, poss, vels = evolve!(p, bt, 10.0)
enableplotting()
plot_billiard(bt)
animate_evolution(p, bt, 200)
