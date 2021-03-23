"""
`kfr(R, A, K)`

Kim-Forger activator-repressor inhibition.
"""
kfr(R, A, K) = (A - R - K + sqrt((A - R - K)^2 + 4*A*K)) / (2*A)


"""
`kfr(R, A)`

Kim-Forger activator-repressor inhibition for K = 0.
"""
kfr(R, A) = R < A ? 1 - R / A : 0.0


"""
`model!(du, u, p, t)`

"""
function model!(du, u, p, t)

    # Species
    M, P, R = u

    # Parameters
    vm, A, dm, vp, dp, vr, dr, I = p

    # Equations
    du[1] = dM = vm*kfr(R, A) - dm*M + I
    du[2] = dP = vp*M - dp*P
    du[3] = dR = vr*P - dr*R

end


"""
`model06!(du, u, p, t)`

"""
function model06!(du, u, p, t)

    # Species
    M, P, R = u

    # Parameters
    A, dm, dp, dr, I = p

    # Equations
    du[1] = dM = kfr(R, A) - dm*M + I
    du[2] = dP = M - dp*P
    du[3] = dR = P - dr*R

end


"""
`model05!(du, u, p, t)`

"""
function model05!(du, u, p, t)

    # Species
    M, P, R = u

    # Parameters
    A, dm, d, I = p

    # Equations
    du[1] = dM = kfr(R, A) - dm*M + I
    du[2] = dP = M - d*P
    du[3] = dR = P - d*R

end


"""
`model04!(du, u, p, t)`

"""
function model04!(du, u, p, t)

    # Species
    M, P, R = u

    # Parameters
    A, d, I = p

    # Equations
    du[1] = dM = kfr(R, A) - d*M + I
    du[2] = dP = M - d*P
    du[3] = dR = P - d*R

end


"""
`loadmodel(name)`

Load a `BCModel` specified by `name`.
"""
function loadmodel(name)
    name = Symbol(name)

    if name == :model
        species = (:M, :P, :R)
        parameters = (:vm, :A, :dm, :vp, :dp, :vr, :dr, :I)
        outfun = sol -> kfr.(sol[3, :], sol.prob.p[2])
        prob = ODEProblem(
            model!,
            [0.1, 0.1, 0.1],
            (0.0, 250.0),
            [10.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0]
        )
        return BCModel(species, parameters, outfun, prob)

    elseif name == :modelsde
        species = (:M, :P, :R)
        parameters = (:vm, :A, :dm, :vp, :dp, :vr, :dr, :I, :noise)
        outfun = sol -> kfr.(sol[3, :], sol.prob.p[2])
        prob = SDEProblem(
            model!,
            (du, u, p, t) -> du[:] .= p[end],
            [0.1, 0.1, 0.1],
            (0.0, 250.0),
            [10.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.1]
        )
        return BCModel(species, parameters, outfun, prob)

    elseif name == :model06
        species = (:M, :P, :R)
        parameters = (:A, :dm, :dp, :dr, :I)
        outfun = sol -> kfr.(sol[3, :], sol.prob.p[1])
        prob = ODEProblem(
            model06!,
            [0.1, 0.1, 0.1],
            (0.0, 250.0),
            [25.0, 0.15, 0.15, 0.15, 0.0]
        )
        return BCModel(species, parameters, outfun, prob)

    elseif name == :model06sde
        species = (:M, :P, :R)
        parameters = (:A, :dm, :dp, :dr, :I, :noise)
        outfun = sol -> kfr.(sol[3, :], sol.prob.p[1])
        prob = SDEProblem(
            model06!,
            (du, u, p, t) -> du[:] .= p[end],
            [0.1, 0.1, 0.1],
            (0.0, 250.0),
            [25.0, 0.15, 0.15, 0.15, 0.0, 0.0]
        )
        return BCModel(species, parameters, outfun, prob)

    elseif name == :model05
        species = (:M, :P, :R)
        parameters = (:A, :dm, :d, :I)
        outfun = sol -> kfr.(sol[3, :], sol.prob.p[1])
        prob = ODEProblem(
            model05!,
            [0.1, 0.1, 0.1],
            (0.0, 250.0),
            [25.0, 0.15, 0.15, 0.0]
        )
        return BCModel(species, parameters, outfun, prob)

    elseif name == :model05sde
        species = (:M, :P, :R)
        parameters = (:A, :dm, :d, :I, :noise)
        outfun = sol -> kfr.(sol[3, :], sol.prob.p[1])
        prob = SDEProblem(
            model05!,
            (du, u, p, t) -> du[:] .= p[end],
            [0.1, 0.1, 0.1],
            (0.0, 250.0),
            [25.0, 0.15, 0.15, 0.0, 0.0]
        )
        return BCModel(species, parameters, outfun, prob)

    elseif name == :model04
        species = (:M, :P, :R)
        parameters = (:A, :d, :I)
        outfun = sol -> kfr.(sol[3, :], sol.prob.p[1])
        prob = ODEProblem(
            model04!,
            [0.1, 0.1, 0.1],
            (0.0, 250.0),
            [25.0, 0.15, 0.0]
        )
        return BCModel(species, parameters, outfun, prob)

    elseif name == :model04sde
        species = (:M, :P, :R)
        parameters = (:A, :d, :I, :noise)
        outfun = sol -> kfr.(sol[3, :], sol.prob.p[1])
        prob = SDEProblem(
            model04!,
            (du, u, p, t) -> du[:] .= p[end],
            [0.1, 0.1, 0.1],
            (0.0, 250.0),
            [25.0, 0.15, 0.0, 0.0]
        )
        return BCModel(species, parameters, outfun, prob)

    else
        throw("Unknown model!")

    end

end