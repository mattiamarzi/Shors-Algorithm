To execute the code you first have to create the following environment:

conda create -n my-quantum-env python=3.8 conda activate my-quantum-env
pip install qiskit
pip install qiskit[visualization]
pip install spyder
pip install numpy==1.19.3

I suggest not to use Python 3.9 and Numpy 1.19.4 because they are not fully compatible.