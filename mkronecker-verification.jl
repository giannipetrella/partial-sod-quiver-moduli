using QuiverTools
using QuiverTools: Chern_character_universal_bundle, Chern_class_line_bundle, Chern_character_from_classes
using QuiverTools: dual_Chern_character as dual



function mKronecker_cohomology_vanishing(M::QuiverModuliSpace;
    chi=QuiverTools.extended_gcd(M.d)[2]
    )

    @info "Checking the following case:\n\n"

    @info "Quiver: $(M.Q)"
    @info "Dimension vector: $(M.d)"
    @info "Stability parameter: $(M.theta)"
    @info "Linearisation: $(chi)\n\n"

    @info "The resulting quiver moduli has:\n\n"

    HH0 = sum(Hodge_diamond(M))
    @info "Dimension: $(QuiverTools.dimension(M))"
    @info "HH_{0}(M) = $HH0\n\n"

    if (nvertices(M.Q) != 2) || (M.d != [2, 3])
        @error("This script is only for m-Kronecker quivers with dimension vector (2, 3). Use `verification.jl` instead.")
    end
    # Chern characters
    cUi = [Chern_character_universal_bundle(M, i; chi=chi) for i in 1:nvertices(M.Q)]
    cUidual = [dual(M, u; chi=chi) for u in cUi]


    Ui_Euler = all(integral(M, v; chi=chi) == 0 for v in cUidual)

    @info "H^{k}(U^∨_i) = 0: $Ui_Euler\n\n"

    return nothing
end

function main()


    a = [2, -1]
    for m in 3:5
        Q = mKronecker_quiver(m)
        M = QuiverModuliSpace(Q, [2, 3])
        mKronecker_cohomology_vanishing(M; chi=a)
    end

end

main()
