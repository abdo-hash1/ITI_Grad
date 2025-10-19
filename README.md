# **QPSK Transceiver Simulation in MATLAB**

## **1\. Overview**

This project is a comprehensive end-to-end simulation of a QPSK (Quadrature Phase-Shift Keying) digital communication system, implemented in MATLAB. It models the entire signal chain, including a transmitter, a channel with realistic impairments, and a sophisticated receiver capable of synchronization and data recovery.  
The primary purpose of this simulation is to analyze the system's Bit Error Rate (BER) performance under various channel conditions, specifically in the presence of Additive White Gaussian Noise (AWGN) and Carrier Frequency Offset (CFO). It serves as a practical tool for understanding and visualizing key concepts in digital signal processing and communications theory.  
[Image of a QPSK constellation diagram](https://encrypted-tbn2.gstatic.com/licensed-image?q=tbn:ANd9GcT6AF6Fra5s4WyDEMLPA0wzEJJgTOW8T27wYF3RWqHNAWy_1lRslsqDu24TYP2a8eBkHWsEc37D9AfHRfCdZ6XZVYiC0KknSrCpqsexRq8bMx2Haq8)

## **2\. System Architecture & Features**

The simulation is built in distinct blocks that model the flow of a signal in a real-world communication system.

### **Key Features Implemented:**

* **QPSK Modulation/Demodulation:** Full implementation of QPSK symbol mapping and demapping.  
* **Pulse Shaping:** Utilizes a Root-Raised Cosine (RRC) filter to shape the transmitted signal, minimizing Inter-Symbol Interference (ISI) at the receiver.  
* **Channel Modeling:**  
  * **Carrier Frequency Offset (CFO):** Simulates frequency mismatch between the transmitter and receiver oscillators.  
  * **AWGN:** Incorporates a correctly implemented AWGN channel based on the desired Energy per Bit to Noise Power Spectral Density ratio (Eb/No).  
* **Advanced Receiver Synchronization:**  
  * **Matched Filtering:** An RRC filter at the receiver, matched to the transmitter's filter, to maximize the signal-to-noise ratio.  
  * **Coarse CFO Estimation:** A preamble-based autocorrelation algorithm to estimate and correct large frequency offsets.  
* **Performance Analysis:**  
  * Calculates the BER of the system across a range of SNR values.  
  * Generates plots comparing the simulated BER against the theoretical QPSK BER curve to validate the model's accuracy.

## **3\. Core Algorithms Explained**

### **a. Pulse Shaping & Matched Filtering**

To prevent spectral overlap and ISI, the baseband signal is shaped using an **RRC filter**. The key property of this approach is that when the signal passes through an identical matched RRC filter at the receiver, the combined response is that of a **Raised Cosine (RC) filter**. This RC filter satisfies the Nyquist criterion for zero ISI, meaning that at the correct sampling instant, the contribution from all other symbols is zero.

### **b. Coarse CFO Estimation**

A significant frequency offset can prevent the receiver from locking onto the signal. This simulation uses a data-aided, feedforward approach to estimate the coarse CFO.

1. A known preamble sequence is prepended to the data packet.  
2. At the receiver, the modulation of the preamble is removed by multiplying it with the complex conjugate of the known sequence.  
3. The resulting signal is a complex sinusoid whose frequency is equal to the CFO.  
4. An autocorrelation is performed on this signal. The phase of the correlation result is directly proportional to the frequency offset, which can then be calculated and corrected.

