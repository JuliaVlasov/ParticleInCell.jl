"""

    Implements split operators for linear pendulum 

    Solve linear pendulum problem: ``\\frac{delta1}{dt} = v ``, 

    ``\\frac{delta2}{dt} = - ω^2 x``. The exact solution is
    ``x(t)=  x(0)   cos(ωt) + \frac{v(0)}{ω}sin(ωt),```
    ``v(t)= -x(0) ω sin(ωt) + v(0)cos(ωt) ``

"""
function linear_pendulum( t_final, x0, v0 )
    x = x0 * cos( t_final ) + v0 * sin( t_final )
    v = -x0 * sin( t_final ) + v0 * cos( t_final )
    x, v
end

function run( nstep )

    x0 = 1.0       # initial x  for order checking
    v0 = 2.0       # initial v  for order checking 
    t_final = 1.0  # final time for order checking 

    x_exact, v_exact = linear_pendulum( t_final, x0, v0 )

    # do iterations with smallest time step
    dt = t_final/nstep

    x = x0
    v = v0

    for istep in 1:nstep
        x += v * dt 
        v -= x * dt
    end
  
    # return  mean square error
    sqrt( (x - x_exact)^2 + (v - v_exact)^2 )


end 


@test run( 10000 ) ≈ 0.0 atol=1e-3
