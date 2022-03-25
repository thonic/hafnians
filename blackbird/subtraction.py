import numpy as np
import strawberryfields as sf
from strawberryfields.ops import *

prog = sf.Program(2)

with prog.context as q:
    Squeezed(1, 0) | q[0]
    BSgate(0.1, 0) | (q[0], q[1])
    MeasureFock(select=1) | q[1]
    Squeezed(1, np.pi) | q[0]

eng = sf.Engine("fock", backend_options={"cutoff_dim": 10})
result = eng.run(prog)
state = result.state
print(state.trace())
print(np.diagonal(state.reduced_dm([0])))

# result:
# [0.64805427+0.j,
# 0.        +0.j,
# 0.18794405+0.j,
# 0.        +0.j,
# 0.08175928+0.j,
# 0.        +0.j,
# 0.03951873+0.j,
# 0.        +0.j,
# 0.02005664+0.j,
# 0.        +0.j]
