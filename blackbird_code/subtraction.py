import numpy as np
import strawberryfields as sf
from strawberryfields.ops import *


def state_generation():
    prog = sf.Program(3)

    with prog.context as q:
        Squeezed(0.8814, 0) | q[0]
       # S2gate(1,0)| (q[0], q[1])
        # S2gate(1.1,0)| (q[0], q[1]
        Dgate(1.4152,np.pi) | q[0]

        BSgate(0.001, 0) | (q[0], q[1])
        MeasureFock(select=1) | q[1]
        Dgate(1.4144,0) | q[0]

        BSgate(0.001, 0) | (q[0], q[2])
        MeasureFock(select=1) | q[2]
       # Dgate(0.4855, -2.3276) | q[0]
        Sgate(0.8814, np.pi) | q[0]
        #MeasureFock(select=1) | q[1]

    # eng = sf.Engine("gaussian")
    # result = eng.run(prog)
    # return result.state_data  # (tuple(mu, cov))

    eng = sf.Engine("fock", backend_options={"cutoff_dim": 25})
    result = eng.run(prog)
    state = result.state
    return np.array([state.trace(), ]), np.diagonal(state.reduced_dm([0]))[:6]