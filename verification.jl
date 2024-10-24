using QuiverTools
using QuiverTools: Chern_character_universal_bundle, Chern_class_line_bundle, Chern_character_from_classes
using QuiverTools: dual_Chern_character as dual

function is_strongly_amply_stable(M::QuiverModuliSpace)
    dest = QuiverTools.all_destabilizing_subdimension_vectors(M.d, M.theta, M.denom)
    return all(e -> Euler_form(M.Q, e, M.d - e) <= -2, dest)
end

function cohomologies_vanishing(M::QuiverModuliSpace;
    chi=QuiverTools.extended_gcd(M.d)[2],
    r=gcd(canonical_stability(M.Q, M.d))
    )

    # check if assumptions hold.
    # if not, throw a warning
    if !(is_acyclic(M.Q) && is_coprime(M.d, M.theta) && is_strongly_amply_stable(M))
        @warn "The assumptions for the verification do not hold.\n\n"
    end


    @info "Checking the following case:\n\n"

    @info "Quiver: $(M.Q)"
    @info "Dimension vector: $(M.d)"
    @info "Stability parameter: $(M.theta)"
    @info "Linearisation: $(chi)\n\n"

    @info "The resulting quiver moduli has:\n\n"

    HH0 = sum(Hodge_diamond(M))
    @info "Dimension: $(QuiverTools.dimension(M))"
    @info "Index: $r"
    @info "HH_{0}(M) = $HH0\n\n"


    # Teleman weights
    HN = all_HN_types(M, unstable=true)
    UidualUj = all_weights_endomorphisms_universal_bundle(M)
    Ui = all_weights_universal_bundle(M; chi=chi)
    Uidual = Dict((hn, -1 .* Ui[hn]) for hn in keys(Ui))
    H = all_weights_irreducible_component_canonical(M)
    H = Dict((hn, -1 .* H[hn]) for hn in keys(H))

    eta = all_Teleman_bounds(M)

    # Chern characters
    cUi = [Chern_character_universal_bundle(M, i; chi=chi) for i in 1:nvertices(M.Q)]
    cUidual = [dual(M, u; chi=chi) for u in cUi]
    cH = Int.(1//r * canonical_stability(M.Q, M.d))
    cH = Chern_character_line_bundle(M, cH; chi=chi)
    cHdual = dual(M, cH; chi=chi)


    QA_predicted = nvertices(M.Q) + 1
    QB_predicted = r * (nvertices(M.Q) + 1)


    @info "Question A predicts $QA_predicted exceptional objects."
    @info "Question B predicts $QB_predicted exceptional objects."
    @info "HH_{0}(M) = $HH0 is the maximum allowed number of exceptional objects.\n\n"

    if QA_predicted > HH0
        @info "Both Questions A and B predict too many exceptional objects.\n\n"

    elseif QB_predicted > HH0
        @info "Question B predicts too many exceptional objects.
      However, Question A might still have a positive answer.\n\n"
    end


    UidualUj_Teleman = map(s -> all(hntype -> maximum(UidualUj[hntype]) - s * H[hntype] < eta[hntype], HN), 1:r - 1)
    UidualUj_Euler = map(r -> all(integral(M, u*v*cHdual^r; chi=chi) == 0 for u in cUidual, v in cUi), 1:r - 1)

    Hdual_Teleman = map(s -> all(hntype -> - s * H[hntype] < eta[hntype], HN), 1:r - 1)
    Hdual_Euler = map(r -> integral(M, cHdual^r; chi=chi) == 0, 1:r - 1)

    Uidual_Teleman = map(s -> all(hntype -> maximum(Uidual[hntype]) - s * H[hntype] < eta[hntype], HN), 0:r)
    Uidual_Euler = map(r -> all(integral(M, u*cHdual^r; chi=chi) == 0 for u in cUidual), 0:r)

    Ui_Teleman = map(s -> all(hntype -> maximum(Ui[hntype]) - s * H[hntype] < eta[hntype], HN), 0:r)
    Ui_Euler = map(r -> all(integral(M, v*cHdual^r; chi=chi) == 0 for v in cUi), 0:r)


    UidualUj_vanishing  = map(s -> (UidualUj_Teleman[s]   && UidualUj_Euler[s] ) || (UidualUj_Teleman[r - s]  && UidualUj_Euler[r - s]),   1:r - 1)
    Hdual_vanishing     = map(s -> (Hdual_Teleman[s]      && Hdual_Euler[s]    ) || (Hdual_Teleman[r - s]     && Hdual_Euler[r - s]),      1:r - 1)
    Uidual_vanishing    = map(s -> (Uidual_Teleman[s + 1] && Uidual_Euler[s + 1])|| (Ui_Teleman[r - s + 1]    && Ui_Euler[r - s + 1]),     0:r - 1)
    Ui_vanishing        = map(s -> (Ui_Teleman[s + 1]     && Ui_Euler[s + 1]   ) || (Uidual_Teleman[r - s + 1]&& Uidual_Euler[r - s + 1]), 0:r - 1)

    for s in 1:r-1
        @info "H^{k}(U^∨_i ⨷ U_j ⨷ O(-$(s)H)) = 0: $(UidualUj_vanishing[s] ? "true" : "?")"
    end
    @info "\n"
    for s in 1:r-1
        @info "              H^{k}(O(-$(s)H)) = 0: $(Hdual_vanishing[s] ? "true" : "?")"
    end
    @info "\n"
        @info "               H^{k}(U^∨_i) = 0: $(Uidual_vanishing[1] ? "true" : "?")"
    for s in 1:r-1
        @info "      H^{k}(U^∨_i ⨷ O(-$(s)H)) = 0: $(Uidual_vanishing[s+1] ? "true" : "?")"
    end
    @info "\n"
    for s in 1:r-1
        @info "        H^{k}(U_j ⨷ O(-$(s)H)) = 0: $(Ui_vanishing[s+1] ? "true" : "?" )"
    end
    print("\n\n")

    QA = Uidual_vanishing[1]
    @info "For the given input, Question A is $(QA)."

    if r >= 2
        QB = all(UidualUj_vanishing) && all(Hdual_vanishing) && all(Uidual_vanishing) && all(Ui_vanishing[2:end])
        @info "For the given input, Question B is $(QB)."
    else
        @info "Question B is equivalent to Question A, as r = 1."
    end

    print("\n\n")
    if QA
        collection = "O, "
        collection *= reduce(*, map(i -> "U_{$(i)}, ", 1:nvertices(M.Q)))

        if r >= 2 && QB
            for s in 1:r-1
                collection *= "O($(s)H), " * reduce(*, map(i -> "U_{$i}($(s)H), ", 1:nvertices(M.Q)))
            end
        end
        # remove the last two characters from the string collection
        collection = collection[1:end-2]

        @info "An exceptional collection is given by $collection.\n\n"
    end

    if !QA
        return nothing
    end
    @info "Now checking if the exceptional collection is strongly exceptional:\n\n"

    if all(hntype -> H[hntype] <= 0, HN)
        @info "All the Harder--Narasimhan weigths of H are negative, so for all s > 0 we have\n
            H^{≥ 1}(U^∨_i ⨷ U_j ⨷ O(sH)) = 0,
            H^{≥ 1}(O(sH))               = 0,
            H^{≥ 1}(U^∨_i ⨷ O(sH))       = 0, and
            H^{≥ 1}(U_j ⨷ O(sH))         = 0.\n\n"

        # H having negative weights leaves with only H^{≥ 1}(U_i) = 0 to check
        last_check = Ui_Teleman[1]
        if last_check
            @info "Moreover, Teleman quantization gives\n
                            H^{≥ 1}(U_i) = 0.\n\n"
            @info "Therefore, the exceptional collection is strongly exceptional.\n\n"
        else
            @info "Howver, Teleman quantization does not give H^{≥ 1}(U_i) = 0.\n\n"
            @info "Therefore, we can't conclude that the exceptional collection is strongly exceptional.\n\n"
        end

        return nothing
    end

    # if the weights of H are not all negative, we must check
    UidualUj_Teleman_op = map(s -> all(hntype -> maximum(UidualUj[hntype]) + s * H[hntype] < eta[hntype], HN), 1:r - 1)
    UidualUj_Euler_op = map(r -> all(integral(M, u*v*cH^r; chi=chi) == 0 for u in cUidual, v in cUi), 1:r - 1)

    Hdual_Teleman_op = map(s -> all(hntype -> + s * H[hntype] < eta[hntype], HN), 1:r - 1)
    Hdual_Euler_op = map(r -> integral(M, cH^r; chi=chi) == 0, 1:r - 1)

    Uidual_Teleman_op = map(s -> all(hntype -> maximum(Uidual[hntype]) + s * H[hntype] < eta[hntype], HN), 0:r)
    Uidual_Euler_op = map(r -> all(integral(M, u*cH^r; chi=chi) == 0 for u in cUidual), 0:r)

    Ui_Teleman_op = map(s -> all(hntype -> maximum(Ui[hntype]) + s * H[hntype] < eta[hntype], HN), 0:r)
    Ui_Euler_op = map(r -> all(integral(M, v*cH^r; chi=chi) == 0 for v in cUi), 0:r)


    # applying Serre duality
    UidualUj_op_vanishing = map(s -> (UidualUj_Teleman_op[s]  && UidualUj_Euler_op[s])   || (UidualUj_Teleman_op[r - s]  && UidualUj_Euler_op[r - s]),  1:r - 1)
    Hdual_op_vanishing    = map(s -> (Hdual_Teleman_op[s]     && Hdual_Euler_op[s])      || (Hdual_Teleman_op[r - s]     && Hdual_Euler_op[r - s]   ),  1:r - 1)
    Uidual_op_vanishing   = map(s -> (Uidual_Teleman_op[s + 1]&& Uidual_Euler_op[s + 1]) || (Ui_Teleman_op[r - s + 1]    && Ui_Euler_op[r - s + 1]),    0:r - 1)
    Ui_op_vanishing       = map(s -> (Ui_Teleman_op[s + 1]    && Ui_Euler_op[s + 1])     || (Uidual_Teleman_op[r - s + 1]&& Uidual_Euler_op[r - s + 1]),0:r - 1)


    for s in 1:r-1
        @info "H^{k}(U^∨_i ⨷ U_j ⨷ O($(s)H)) = 0: $(UidualUj_op_vanishing[s] ? "true" : "?")"
    end
    for s in 1:r-1
        @info "               H^{k}(O($(s)H)) = 0: $(Hdual_op_vanishing[s] ? "true" : "?")"
    end
    for s in 1:r-1
        @info "      H^{k}(U^∨_i ⨷ O($(s)H)) = 0: $(Uidual_op_vanishing[s] ? "true" : "?")"
    end
        @info "                   H^{k}(U_j) = 0: $(Ui_op_vanishing[1] ? "true" : "?")"
    for s in 1:r-1
        @info "        H^{k}(U_j ⨷ O($(s)H)) = 0: $(Ui_op_vanishing[s + 1] ? "true" : "?")"
    end

    if r == 1
        strong = Ui_op_vanishing[1]
    else
        strong = all(UidualUj_op_vanishing) && all(Hdual_op_vanishing) && all(Uidual_op_vanishing[2:end]) && all(Ui_op_vanishing)
    end

    @info "Strongly exceptional: $(strong)"

    return nothing
end





function main()


    # del Pezzo surfaces

    # P^2
    #
    Q = mKronecker_quiver(3);
    d = [1, 1];
    a = [2, -1];
    M = QuiverModuliSpace(Q, d);
    cohomologies_vanishing(M; chi=a)

    # P^1 x P^1
    #
    Q = Quiver("1--3,2--3");
    d = [1, 1, 1];
    a = [1, 1, -1];
    M = QuiverModuliSpace(Q, d);
    cohomologies_vanishing(M; chi=a)

    # P^2 blown up at 1 point
    #
    Q = Quiver("1-2,1-3,2--3");
    d = [1, 1, 1];
    a = [1, 1, -1];
    M = QuiverModuliSpace(Q, d);
    cohomologies_vanishing(M; chi=a)

    # P^2 blown up at 2 points
    #
    Q = Quiver("1-3,1-4,2-3,2-4,3-4");
    d = [1, 1, 1, 1];
    a = [1, 1, 0, -1];
    M = QuiverModuliSpace(Q, d);
    cohomologies_vanishing(M; chi=a)

    # P^2 blown up at 3 points
    #
    Q = Quiver("1-4,2-4,3-4,1-5,2-5,3-5");
    d = [1, 1, 1, 1, 1];
    a = [1, 1, 1, -1, -1];
    M = QuiverModuliSpace(Q, d);
    cohomologies_vanishing(M; chi=a)

    # P^2 blown up at 4 points
    #
    Q = Quiver("1-6,2-6,3-6,4-6,5-6");
    d = [1, 1, 1, 1, 1, 2];
    a = [1, 1, 1, 1, 1, -2];
    M = QuiverModuliSpace(Q, d);
    cohomologies_vanishing(M; chi=a)




    # Framed quiver for Fano 3-fold
    #
    Q = Quiver("1-2,1-3,2---3");
    d = [1, 1, 1];
    a = [1, 1, -1];
    M = QuiverModuliSpace(Q, d);
    cohomologies_vanishing(M; chi=a)


    # Fano 5-fold
    #
    Q = Quiver("1-2,1--3,2---3");
    d = [1, 1, 1];
    a = [0, 2, -1];
    M = QuiverModuliSpace(Q, d);
    cohomologies_vanishing(M; chi=a)

    # Recurring example in the paper

    Q = mKronecker_quiver(3);
    d = [3, 4];
    a = [3, -2];
    M = QuiverModuliSpace(Q, d);
    cohomologies_vanishing(M; chi=a)


    # m-Kronecker quiver for various values of m

    m_check = 5
    d = [2, 3];
    a = [2, -1]
    @warn "this is extremely slow!"
    for m in 3:m_check
       Q = mKronecker_quiver(m);
       M = QuiverModuliSpace(Q, d)
       cohomologies_vanishing(M; chi=a);
    end

end


main()
