# -*- coding: utf-8 -*-
"""
Shor's Algorithm

@author: mattia
"""

import numpy as np
from qiskit import QuantumCircuit, Aer, execute
from qiskit.visualization import plot_histogram
import pandas as pd
from fractions import Fraction
import time

#Give a starting time to measure execution time
start = time.time()

#Define the function c_7mod15 which returns controlled-U gate for m repeated x times
#c_7mod15 will be a 4 qubit unitary matrix controlled by a 5th qubit which will be added
def c_7mod15(power):
    "Controlled multiplication by 7 mod 15"
    U = QuantumCircuit(4) #crea un circuito a 4 qubits
    for iteration in range(power):
        U.swap(2,3)
        U.swap(1,2)
        U.swap(0,1)
        for q in range(4):
            U.x(q) #Apply X gate to the q-th qubit
    U = U.to_gate() #Define a custom gate
    U.name = "%i^%i mod 15" %(7, power) #Give a name to the custom gate
    c_U = U.control() #Make the custom gate controlled
    return c_U


#Define the circuit for the inverse QFT
def qft_dagger(n):
    "inverse quantum Fourier transform to the first n qubits"
    qc = QuantumCircuit(n)
    for qubit in range(n//2): #The swaps are done at the beginning because this is the inverse QFT; if n is odd then the central term doesn't need to be swapped
        qc.swap(qubit, n-qubit-1) #n-qubit indicates the (n-qubit)-th qubit, note that we start counting the qubits from 0 to n-1
    for j in range(n):
        for m in range(j):
            qc.cu1(-2*np.pi/float(2**(j-m+1)), m, j) #the function qc.cu1(x, q[0], q[1]) applies the rotation (which is diagonal with all 1 except for the last term which is e^(ix))  to the target q[1], only if the control q[0] is 1
        qc.h(j)
    qc.name = "QFT^(-1)"
    return qc
qft_dagger(4).draw(output='mpl', filename='Inverse QFT')


#Let's use 8 counting qubits (more counting qubits = bigger phase predictable) plus 4 qubits for U to act on
n_count = 8
qc = QuantumCircuit(n_count + 4, n_count) #The second argument indicates the number of bits (which will be used to store information of the measurements)
for q in range(n_count): #Applying the Hadamard gate to the counting qubits (at the begininning they are initialized to 0)
    qc.h(q)
qc.x(3+n_count) #The last ancilla qubit becomes 1, in this way the ancilla register is initialized to 0001
for q in range(n_count):
    qc.append(c_7mod15(2**(q)), [q] + [i+n_count for i in range(4)]) #Applies the gate c_7mod15(2**q) to the ancilla qubits [1+n_count], [2+n_count], [3+n_count], [4+n_count] utilizing the qubit [q] as control
qc.append(qft_dagger(n_count), range(n_count)) #Applies the gate qft_dagger(n_count) to the n_count qubits
qc.measure(range(n_count), range(n_count)) #Measure the n_count qubits storing information to the n_count bits
qc.draw(output='mpl', filename='Shor circuit implementation') #Displays the image of the circuit saving it as "Shor's circuit implementation"


#Let's see what results we measure
backend = Aer.get_backend('qasm_simulator') #Connecting to the Aer simulator of IBM
results = execute(qc, backend, shots = 2048).result() #Execute the circuit qc using the simulator backend 2048 time
counts = results.get_counts() #Creates an array: "results: number of times I got that result"
fig = plot_histogram(counts)
fig.savefig('istogramma.png', bbox_inches = 'tight') #The bbox command is necessary to make the figure of the right dimension
rows, measured_phases = [], []
for output in counts:
    decimal = int(output, 2) #convert base 2 string to decimal
    phase = decimal/(2**n_count)
    measured_phases.append(phase)
    rows.append(["%s(bin) = %i(dec)" %(output, decimal),
                 "%i/%i = %.2f" %(decimal, 2**n_count, phase)])

#Print the rows in a table 
headers = ["Register Output", "Phase"]
df = pd.DataFrame(rows, columns = headers)
print(df)

#Creating a table with the guesses for P
rows = []
for phase in measured_phases:
    frac =  Fraction(phase).limit_denominator(15)
    rows.append([phase, "%i/%i" %(frac.numerator, frac.denominator), frac.denominator])
headers = ["Phase", "Fraction", "Guesses for P"]
df = pd.DataFrame(rows, columns = headers)
print(df)

end = time.time()
print("Execution time = %5.3f" %(end - start))
