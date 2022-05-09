import numpy as np
import strawberryfields as sf
from strawberryfields.ops import *


def state_generation():
    prog = sf.Program(3)

    with prog.context as q:
        Squeezed(1, 0) | q[0]
       # S2gate(1,0)| (q[0], q[1])
        # S2gate(1.1,0)| (q[0], q[1]
        Dgate(1.4144,0) | q[0]
        BSgate(0.001, 0) | (q[0], q[1])
        MeasureFock(select=1) | q[1]
        Dgate(1.4144,np.pi) | q[0]
        BSgate(0.001, 0) | (q[0], q[2])
        MeasureFock(select=1) | q[2]
        Sgate(1, np.pi) | q[0]
     #   MeasureFock(select=1) | q[1]


    eng = sf.Engine("fock", backend_options={"cutoff_dim": 50})
    result = eng.run(prog)
    state = result.state
    return np.array([state.trace(), ]), np.diagonal(state.reduced_dm([0]))[:6]

    # # result:
    # 0.9999997939494688
    # [0.64805427+0.j 0.        +0.j 0.18794405+0.j 0.        +0.j]
