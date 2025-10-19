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
  * **Fine Carrier Tracking:** A decision-directed **Phase-Locked Loop (PLL)** (Costas Loop) to track and correct any residual frequency and phase offsets after the coarse correction.  
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

### **c. Phase-Locked Loop (PLL) for Fine Tracking**

After coarse correction, a **Costas Loop PLL** is used to track small, residual frequency errors and the signal's phase.

* **Operation:** The loop operates on a symbol-by-symbol basis. For each incoming symbol, it calculates a phase error by comparing the symbol's phase against the phase of its closest constellation point (the "decision").  
* **Loop Filter:** This error signal is fed into a Proportional-Integral (PI) loop filter, which smooths the error and controls a Numerically Controlled Oscillator (NCO).  
* **NCO:** The NCO generates the phase correction that is applied to the next incoming symbol, effectively "locking" the receiver's phase to the signal's phase.  
* **Data-Aided Handoff:** The loop uses the known preamble symbols for initial fast and accurate locking (data-aided mode) before switching to using its own decisions for tracking the payload data (decision-directed mode).

## **4\. How to Run the Simulation**

1. Ensure all script files (.m) are in the same MATLAB directory or path. This includes the main script and the helper functions (designPulseShapingFilter.m, add\_awgn\_ebno.m, etc.).  
2. Open the main simulation script.  
3. Adjust the parameters in the "**System Parameters**" section at the top of the file:  
   * SNR\_dB\_range: The range of Eb/No values to simulate.  
   * random\_cfo\_value: The actual CFO (in Hz) to introduce in the channel.  
   * num\_of\_packets: The number of packets to simulate for each SNR point. **Increase this value for smoother, more accurate BER curves (e.g., 100 or more).**  
   * Kp, Ki: The proportional and integral gains for the PLL. These may need tuning for different channel conditions.  
4. Run the script.  
5. The simulation progress will be printed to the MATLAB command window, and the final BER plot will be generated in a figure window.

## **5\. Sample Simulation Result**

The following plot shows the output of the simulation, comparing the measured BER against the theoretical performance of QPSK. The close match validates the correctness of the entire simulated transceiver chain, including the synchronization algorithms.  
*The jaggedness of the simulated curve at lower BER values indicates that more packets should be simulated to achieve statistical convergence.*

## **6\. Dependencies**

* MATLAB (tested on version R2020a and later).  
* No special toolboxes are required; all functions use standard MATLAB capabilities.  
* All required helper functions are included in the project files.

*This README was generated to document the QPSK simulation project. It details the system architecture, algorithms, and usage instructions.*