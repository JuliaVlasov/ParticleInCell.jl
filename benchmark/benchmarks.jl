using BenchmarkTools

const suite = BenchmarkGroup()

suite["pendulum"] = BenchmarkGroup(["integration", "runge-kutta"])
n = 10000
simulation = Simulation(n)
suite["pendulum"]["step"] = @benchmarkable step(simulation)
