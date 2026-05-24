from qiskit import QuantumCircuit, transpile
from qiskit_aer import AerSimulator
from qiskit.visualization import plot_histogram
import matplotlib.pyplot as plt
import time
import numpy as np

def AliceBob():
    q = QuantumCircuit(2, 2)
    q.h(0)
    q.cx(0, 1)
    message = "01"
    if message == "00":
        pass
    elif message == "01":
        q.x(0)
    elif message == "10":
        q.z(0)
    elif message == "11":
        q.z(0)
        q.x(0)

    q.cx(0, 1)
    q.h(0)
    q.measure([0, 1], [0, 1])
    return q

def add(qc: QuantumCircuit, a: int, st: int, fin: int) -> None:
    if st > fin or fin >= qc.num_qubits:
        raise ValueError(f"Некорректные индексы: st={st}, fin={fin}, num_qubits={qc.num_qubits}")  
    n = fin - st + 1
    zn = 1 << n
    two_pi = 2 * np.pi
    for i in range(st, fin + 1):
        bit_pos = i - st
        theta = two_pi * a * (1 << bit_pos) / zn
        qc.p(theta, i)

def QFT_once(qc: QuantumCircuit, st: int, fin: int) -> None:
    if st > fin or fin >= qc.num_qubits:
        raise ValueError(f"Некорректные индексы: st={st}, fin={fin}, num_qubits={qc.num_qubits}")
    qc.h(st)
    s = 2
    for i in range(st + 1, fin + 1):
        qc.cp(np.pi / s, i, st)
        s = s << 1

def QFT(qc: QuantumCircuit, st: int = None, fin: int = None) -> None:
    if st is None:
        st = 0
    if fin is None:
        fin = qc.num_qubits - 1
    for i in range(st, fin + 1):
        QFT_once(qc, i, fin)

def IQFT_once(qc: QuantumCircuit, st: int, fin: int) -> None:
    if st > fin or fin >= qc.num_qubits:
        raise ValueError(f"Некорректные индексы: st={st}, fin={fin}, num_qubits={qc.num_qubits}")
    s = 1 << (fin - st)
    for i in range(fin, st, -1):
        qc.cp(-np.pi / s, i, st)
        s = s >> 1
    qc.h(st)

def IQFT(qc: QuantumCircuit, st: int = None, fin: int = None) -> None:
    if st is None:
        st = 0
    if fin is None:
        fin = qc.num_qubits - 1
    for i in range(fin, st - 1, -1):
        IQFT_once(qc, i, fin)

for i in range(2,22,2):
    q = QuantumCircuit(i, i)
    QFT(q, 0, i-1)
    add(q, 500000, 0, i-1)
    IQFT(q, 0, i-1)
    q.measure(range(i), range(i))
    simulator = AerSimulator(method = 'statevector', max_parallel_threads = 1)
    start = time.perf_counter()
    compiled_circuit = transpile(q, simulator)
    result = simulator.run(compiled_circuit, shots=1000).result()
    counts = result.get_counts()
    end = time.perf_counter()
    print(f"{i}: {end - start:.6f}")
