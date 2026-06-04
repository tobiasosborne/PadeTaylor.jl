# Measure VS.1.2's exact per-component errors under the ADR-0028 dispatch, to
# re-baseline the test tolerances on data (not guesses).  Harmonic [cos,-sin],
# order 30, h 0.7 — the VS.1.2 setup.
using PadeTaylor
using PadeTaylor: VectorPadeStepperState, vector_pade_step!
harm(z, y) = [y[2], -y[1]]

st = VectorPadeStepperState{Float64}(0.0, [1.0, 0.0])
vector_pade_step!(st, harm, 30, 0.7)
println("F64  cos err = ", abs(st.y[1] - cos(0.7)), "   -sin err = ", abs(st.y[2] + sin(0.7)))

setprecision(BigFloat, 256) do
    hb = BigFloat("0.7")
    stb = VectorPadeStepperState{BigFloat}(BigFloat(0), BigFloat[1, 0])
    vector_pade_step!(stb, harm, 30, hb)
    println("BF   cos err = ", Float64(abs(stb.y[1] - cos(hb))),
            "   -sin err = ", Float64(abs(stb.y[2] + sin(hb))))
end
