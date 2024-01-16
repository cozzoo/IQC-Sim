import numpy as np
import cmath
import sys
import io
from PyQt5.QtGui import QPixmap
from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QWidget, QPushButton, QTextEdit, QHBoxLayout, QLabel, QComboBox, QMessageBox
from qiskit import QuantumCircuit, transpile, assemble, Aer, execute
from qiskit.visualization import plot_histogram
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvasAgg
from qiskit import IBMQ



class State:
    def __init__(self, n_qubits):
        # TODO task 1.1
        self.n_qubits = n_qubits
        self.state = np.zeros(2**n_qubits)
        self.state[0]=1 #state 0 - first qubit is 0

    def initialize_state(self, qubit_values):
        assert len(qubit_values) == self.n_qubits

        # Convert qubit values to the corresponding index in the state vector
        index = 0
        for i, qubit in enumerate(qubit_values):
            index += qubit * (2 ** (len(qubit_values) - i - 1))

        # Initialize the state vector with all zeros
        self.state = np.zeros(2 ** self.n_qubits)

        # Set the value at the calculated index to 1
        self.state[index] = 1


    def apply_gate(self, gate, n_gate, starting_qubit):
        # TODO task 1.3
         # Initialize the identity matrix for the entire system
        identity_matrix = np.eye(2)
        no_is_before = starting_qubit
        no_is_after = self.n_qubits - starting_qubit - n_gate

        # Calculate the position to apply the gate_matrix
        matrix_i_before = np.eye(2**no_is_before)
        matrix_i_after = np.eye(2**no_is_after)

        unitary_matrix = np.kron(matrix_i_before, gate)
        unitary_matrix = np.kron(unitary_matrix, matrix_i_after)
        self.state = np.dot(self.state, unitary_matrix)




    # Construct the unitary matrix with the gate applied


    def apply_H_gate(self, target_qubit):
        Hadamard = np.sqrt(1/2)*np.matrix('1,1;1,-1')
        self.apply_gate(Hadamard, 1, target_qubit)
        return self.state
        pass
    
    def apply_Y_gate(self, target_qubit):
        Y_gate = np.array([[0, -1j], [1j, 0]])
        self.apply_gate(Y_gate, 1, target_qubit)
        return self.state
        pass
    
    def apply_Z_gate(self, target_qubit):
        Z_gate = np.array([[1, 0],
                   [0, -1]])
        self.apply_gate(Z_gate, 1, target_qubit)
        return self.state
        pass
    
    def apply_T_gate(self, target_qubit):
        T_gate = np.array([[1.        +0.j        , 0.        +0.j        ], [0.        +0.j        , 0.70710678+0.70710678j]])
        self.apply_gate(T_gate, 1, target_qubit)
        return self.state
        pass

    def apply_X_gate(self, target_qubit):
        XMatrix = np.matrix('0,1;1,0')
        self.apply_gate(XMatrix, 1, target_qubit)
        return self.state
        pass

    def apply_S_gate(self, target_qubit):
        SMatrix = np.array([[1, 0], [0, 1j]])
        self.apply_gate(SMatrix, 1, target_qubit)
        return self.state
        pass

    def apply_CNOT_gate(self, target_qubit):
        CNOT_matrix=np.matrix('1,0,0,0;0,1,0,0;0,0,0,1;0,0,1,0')
        self.apply_gate(CNOT_matrix, 2, target_qubit)
        return self.state
        pass

    def apply_CH_gate(self, target_qubit):
        CH_matrix = np.array([[1, 0, 0, 0],
                                     [0, 1/np.sqrt(2),0, 1/np.sqrt(2)],
                                             [0, 1, 0, 0],
                                     [0, 1/np.sqrt(2),0, -1/np.sqrt(2)]])
        self.apply_gate(CH_matrix, 2, target_qubit)
        return self.state
        pass

    def apply_SWAP_gate(self, target_qubit):
        SWAP_matrix=np.array([[1,0,0,0],
                      [0,0,1,0],
                      [0,1,0,0],
                      [0,0,0,1]])
        self.apply_gate(SWAP_matrix, 2, target_qubit)
        return self.state
        pass

    def apply_CY_gate(self, target_qubit):
        CY_matrix=np.array([[ 1.+0.j,  0.+0.j,  0.+0.j,  0.+0.j],[ 0.+0.j,  1.+0.j,  0.+0.j,  0.+0.j],[ 0.+0.j,  0.+0.j,  0.+0.j, -0.-1.j],[ 0.+0.j,  0.+0.j,  0.+1.j,  0.+0.j]])
        self.apply_gate(CY_matrix, 2, target_qubit)
        return self.state
        pass

    def apply_CZ_gate(self, target_qubit):
        CZ_matrix=np.array([[ 1,  0,  0,  0],[ 0,  1,  0,  0],[ 0,  0,  1,  0],[ 0,  0,  0, -1]])
        self.apply_gate(CZ_matrix, 2, target_qubit)
        return self.state
        pass
    
    def apply_CT_gate(self, target_qubit):
        # TODO task 1.5
        CT_matrix=np.array([[1.        +0.j        , 0.        +0.j        , 0.        +0.j        , 0.        +0.j        ],
                             [0.        +0.j        , 1.        +0.j        , 0.        +0.j        , 0.        +0.j        ],
                             [0.        +0.j        , 0.        +0.j        , 1.        +0.j        , 0.        +0.j        ],
                             [0.        +0.j        , 0.        +0.j        , 0.        +0.j        , 0.70710678+0.70710678j]])

        self.apply_gate(CT_matrix, 2, target_qubit)
        return self.state
        pass
    
    
    def apply_CS_gate(self, target_qubit):
        # TODO task 1.5
        CS_matrix=np.array([[1.+0.j, 0.+0.j, 0.+0.j, 0.+0.j],[0.+0.j, 1.+0.j, 0.+0.j, 0.+0.j],[0.+0.j, 0.+0.j, 1.+0.j, 0.+0.j],[0.+0.j, 0.+0.j, 0.+0.j, 0.+1.j]])

        self.apply_gate(CS_matrix, 2, target_qubit)
        return self.state
        pass
    
    def apply_TOFFOLI_gate(self, target_qubit):
        Toffoli_gate = np.array([[1, 0, 0, 0, 0, 0, 0, 0],
                         [0, 1, 0, 0, 0, 0, 0, 0],
                         [0, 0, 1, 0, 0, 0, 0, 0],
                         [0, 0, 0, 1, 0, 0, 0, 0],
                         [0, 0, 0, 0, 1, 0, 0, 0],
                         [0, 0, 0, 0, 0, 1, 0, 0],
                         [0, 0, 0, 0, 0, 0, 0, 1],
                         [0, 0, 0, 0, 0, 0, 1, 0]], dtype=complex)
        self.apply_gate(Toffoli_gate, 3, target_qubit)
        return self.state
        pass
    
          
    def quantum_fourier_transform(self, n):
        """ 
        Create a Quantum Fourier Transform (QFT) matrix for n qubits.
        """
        N = 2**n  # Total number of states
        qft_matrix = np.zeros((N, N), dtype=complex)

        for x in range(N):
            for y in range(N):
                qft_matrix[x, y] = cmath.exp(2j * np.pi * x * y / N) / np.sqrt(N)

        return qft_matrix
    
    
    def apply_QFT_2_gate(self, target_qubit):
        # TODO task 1.5
        QFT_2_matrix=np.array(self.quantum_fourier_transform(2))
        self.apply_gate(QFT_2_matrix, 2, target_qubit)
        return self.state
        pass
    
    

    def produce_measurement(self):

        # Calculate the probabilities of measuring each basis state
        probabilities = np.abs(self.state)**2
        # Generate a random measurement outcome based on probabilities
        measurement = np.random.choice(len(slef.state), p=probabilities)
        # Convert the measurement outcome to a binary string
        measurement_binary = format(measurement, f'0{int(np.log2(len(self.state)))}b')
        # Convert the binary string to a list of qubit values
        measurements = [int(bit) for bit in measurement_binary]

        return measurements
    


class QuantumSimulator(QMainWindow):
    def __init__(self):
        super().__init__()
        self.initUI()
        
        # Initialize the qubit count
        self.qubitCount = 0
        # Data structures to hold qubits and gates
        self.qubits = []
        self.gateSelections = {}  # Dictionary to store gate selections
        
        self.quantumState = np.array([])
    

    def initUI(self):
        # Main Window Setup
        self.setWindowTitle("Quantum Simulator")
        self.centralWidget = QWidget()
        self.setCentralWidget(self.centralWidget)
        self.layout = QVBoxLayout(self.centralWidget)
        
        

        # Add Qubits Section
        self.addQubitBtn = QPushButton('Add Qubit', self)
        self.addQubitBtn.clicked.connect(self.addQubit)
        self.layout.addWidget(self.addQubitBtn)

        # Initial State Section
        self.initialStateTextBox = QTextEdit(self)
        self.initialStateTextBox.setPlaceholderText("Enter initial state here...")
        self.layout.addWidget(self.initialStateTextBox)
        self.setInitialStateBtn = QPushButton('Set Initial State', self)
        self.setInitialStateBtn.clicked.connect(self.setInitialState)
        self.layout.addWidget(self.setInitialStateBtn)

        # Qubit Display Area
        self.qubitLayout = QVBoxLayout()
        self.layout.addLayout(self.qubitLayout)

        # Calculate State Button
        self.calculateStateBtn = QPushButton('Calculate State', self)
        self.calculateStateBtn.clicked.connect(self.calculateState)
        self.layout.addWidget(self.calculateStateBtn)

        # Output Area
        self.outputArea = QTextEdit(self)
        self.outputArea.setReadOnly(True)
        self.layout.addWidget(self.outputArea)
        
        #Run on IBM Hardware
        #self.runIBMbtn = QPushButton('Run on IBM Hardware', self)
        #self.runIBMbtn.clicked.connect(self.runIBM)

        
       

    def addQubit(self):
        
        self.qubitCount += 1
        qubitRow = QHBoxLayout()
        qubitLabel = QLabel(f'Qubit {self.qubitCount}')
        qubitRow.addWidget(qubitLabel)

        qubitGates = []
        for position in range(5):  # Assuming 5 slots per qubit
            gateComboBox = QComboBox()
            gateComboBox.addItems(["-", "Gate X", "Gate Y", "Gate Z", "Gate T", "Gate H", "Gate CNOT", "Gate CH", "Gate SWAP", "Gate S", "Gate CY", "Gate CZ", "Gate CT", "Gate CS","Gate Toffoli", "Gate QFT_2"])  # Example gates
            qubitRow.addWidget(gateComboBox)
            qubitGates.append((position, gateComboBox))

        self.gateSelections[self.qubitCount] = qubitGates
        self.qubits.append(qubitRow)
        self.qubitLayout.addLayout(qubitRow)

    def setInitialState(self):
        # Get initial state from the text box as a string
        initialStateStr = self.initialStateTextBox.toPlainText()

        # Convert the string to a list of integers
        try:
            initialStateList = [int(bit) for bit in initialStateStr if bit in '01']

            # Check if the length of initialStateList matches the qubitCount
            if len(initialStateList) != self.qubitCount:
                self.outputArea.setText("Error: The number of bits in the initial state must match the number of qubits.")
                return

            # Convert the list to a numpy array and store it
            auxState = State(self.qubitCount)
            print(auxState.state)
            self.quantumState = np.array(initialStateList)
            print(self.quantumState)
            auxState.initialize_state(self.quantumState)
            self.quantumState = auxState.state
            print('Initial quantum state is:')
            print(self.quantumState)
            self.outputArea.setText(f"Initial state set to: {self.quantumState}")
        except ValueError:
            self.outputArea.setText("Error: Initial state must be a binary string (containing only 0s and 1s).")
            
    def density_matrix(state):
        dm = np.outer(state, np.conj(state))
        return dm

    def calculateState(self):
        
        availableDiagramQiskit = 1;
        
        #Declare QuantumCircuit from qiskit lib
        self.qc = QuantumCircuit(self.qubitCount)
        
        # TODO: Apply gates to qubits and compute final state
        finalState = "Final State"  # Placeholder
        self.outputArea.setText(finalState)
        
        stateClass = State(self.qubitCount)
        stateClass.state = self.quantumState
            
        #for qubitNum, gates in self.gateSelections.items():
                #for gate in gates:
                    #position, gateComboBox = gate
                    #print(f"Position: {position}, Type of Position: {type(position)}")
                    #print(f"GateComboBox: {gateComboBox}, Type of GateComboBox: {type(gateComboBox)}")
                    #print(gateComboBox.currentText())
                    #selectedGate = gateComboBox.currentText()
                    #print(self.quantumState)
                # Apply the gate function based on the selectedGate
                # Example: if selectedGate == "Gate 1": apply_gate_1_to_qubit(qubitNum, position)
                # You need to define how to apply these
                    # After applying all gates, compute the final state
                # Example: finalState = self.computeFinalState()
                
                # self.outputArea.setText(f"Final state: {finalState}")

            # Define additional methods as needed, such as methods for applying gates and computing the final state

        for position in range(5):
                for qubitNum, gates in self.gateSelections.items():
                    gate = gates[position]
                    gateComboBox = gate[1]
                    selectedGate = gateComboBox.currentText()
                    print(stateClass.state)
                    if selectedGate == "Gate X":
                        stateClass.apply_X_gate(qubitNum-1)
                        print(stateClass.state)
                        self.qc.x(qubitNum-1)
                    if selectedGate == "Gate Y":
                        stateClass.apply_Y_gate(qubitNum-1)
                        print(stateClass.state)
                        self.qc.y(qubitNum-1)                        
                    if selectedGate == "Gate Z":
                        stateClass.apply_Z_gate(qubitNum-1)
                        print(stateClass.state)
                        self.qc.z(qubitNum-1)
                    if selectedGate == "Gate T":
                        stateClass.apply_T_gate(qubitNum-1)
                        print(stateClass.state)
                        self.qc.t(qubitNum-1)
                    if selectedGate == "Gate H":
                            stateClass.apply_H_gate(qubitNum-1)
                            print(stateClass.state)
                            self.qc.h(qubitNum-1)
                    if selectedGate == "Gate CNOT":
                            stateClass.apply_CNOT_gate(qubitNum-1)
                            print(stateClass.state)
                            self.qc.cx(qubitNum-1, qubitNum)
                    if selectedGate == "Gate S":
                            stateClass.apply_S_gate(qubitNum-1)
                            print(stateClass.state)
                            self.qc.s(qubitNum-1)
                    if selectedGate == "Gate CH":
                            stateClass.apply_CH_gate(qubitNum-1)
                            print(stateClass.state)
                            self.qc.ch(qubitNum-1, qubitNum)
                    if selectedGate == "Gate SWAP":
                            stateClass.apply_SWAP_gate(qubitNum-1)
                            print(stateClass.state)
                            self.qc.cs(qubitNum-1, qubitNum)
                    if selectedGate == "Gate CY":
                            stateClass.apply_CY_gate(qubitNum-1)
                            print(stateClass.state)
                            self.qc.cy(qubitNum-1, qubitNum)
                    if selectedGate == "Gate CZ":
                            stateClass.apply_CZ_gate(qubitNum-1)
                            print(stateClass.state)
                            availableDiagramQiskit=0
                    if selectedGate == "Gate CT":
                            stateClass.apply_CT_gate(qubitNum-1)
                            print(stateClass.state)
                            availableDiagramQiskit=0
                    if selectedGate == "Gate CS":
                            stateClass.apply_CS_gate(qubitNum-1)
                            print(stateClass.state)
                            availableDiagramQiskit=0
                    if selectedGate == "Gate QFT_2":
                            stateClass.apply_QFT_2_gate(qubitNum-1)
                            print(stateClass.state)
                            availableDiagramQiskit=0
                    if selectedGate == "Gate Toffoli":
                            stateClass.apply_TOFFOLI_gate(qubitNum-1)
                            print(stateClass.state)
                            availableDiagramQiskit=0
        finalState = np.array2string(stateClass.state)
        self.outputArea.setText(finalState)
        
        if availableDiagramQiskit == 1:
            self.qc.measure_all()
            # Simulate the quantum circuit
            simulator = Aer.get_backend('qasm_simulator')
            result = execute(self.qc, simulator, shots=1024).result()
            counts = result.get_counts()

            # Generate a histogram
            figure = plot_histogram(counts)
            canvas = FigureCanvasAgg(figure)
            canvas.draw()

            # Save the histogram to a buffer
            buf = io.BytesIO()
            canvas.print_png(buf)


            # Display the histogram in a message box
            pixmap = QPixmap()
            pixmap.loadFromData(buf.getvalue(), format="PNG")
            simulation_result = QLabel()
            simulation_result.setPixmap(pixmap)
            buf.close()

            # Create a message box to display the result
            msg_box = QMessageBox()
            msg_box.setWindowTitle("Simulation Result")
            msg_box.layout().addWidget(simulation_result)
            msg_box.exec_()
        
        
    #def runIBM(self):
     #   IBMQ.save_account('9eb467b26d6cf88c9a1ca26414ad5c9228e046859357af45041d31179b135ac83904e8ade6c72fe9691ea08597f9cee53b7fc2dc216674a4bf257dc4880d8e08')
      #  IBMQ.load_account()

       # # Connect to one of the IBM quantum devices
        #provider = IBMQ.get_provider(hub='ibm-q')
        #backend = provider.get_backend('aer_simulator')  # Replace with an available backend

        # Run the same circuit as before
        #qobj = assemble(qc_compiled, backend)
        #job = backend.run(qobj)
        #result = job.result()

        # Get the results and plot
        #counts = result.get_counts(qc)
        #plot_histogram(counts)
        
                
            

def main():
    app = QApplication(sys.argv)
    ex = QuantumSimulator()
    ex.show()
    sys.exit(app.exec_())
       
if __name__ == '__main__':
        main()