export UA

mutable struct UA
  
    ntau :: Int64
    ε    :: Float64
    tau  :: Vector{Float64}
    ltau :: Vector{Float64}
    ftau :: Array{ComplexF64,1}
    ptau :: FFTW.cFFTWPlan{ComplexF64,-1,false,1}
    pl   :: Array{ComplexF64,2}
    ql   :: Array{ComplexF64,2}
    r    :: Array{ComplexF64,2}
    r̃    :: Array{ComplexF64,2}
    rtau :: FFTW.cFFTWPlan{ComplexF64,-1,false,2}

    function UA( ntau, ε, nbpart )

        dtau = 2π / ntau
        
        ltau  = zeros(Float64, ntau)
        ltau .= vcat(0:ntau÷2-1, -ntau÷2:-1) 
        
        tau   = zeros(Float64, ntau)
        tau  .= [ i*dtau for i=0:ntau-1 ]

        ftau  = zeros(ComplexF64, ntau)

        ptau  = FFTW.plan_fft(ftau)

        pl = zeros(ComplexF64, (ntau, nbpart))
        ql = zeros(ComplexF64, (ntau, nbpart))

        r  = zeros(ComplexF64, (ntau,2))
        r̃  = zeros(ComplexF64, (ntau,2))

        rtau = plan_fft(r,1)

        new( ntau, ε, tau, ltau, ftau, ptau, pl, ql, r, r̃, rtau )

    end

end
