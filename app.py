from flask import Flask, render_template, url_for, request, redirect

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('agg')
import scipy.linalg as LA
import cmath
import math
import os
from werkzeug.utils import secure_filename

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = 'uploads'
os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)

def MEMO_ESPRIT(Voltage,fs):
#     inputSignal = pd.read_csv("abc.csv",header=None) #reading uploaded input csv file
#     inputDataSamples=inputSignal.iloc[0,:]  #array containing inputsignal)
    inputDataSamples=np.transpose(Voltage)
    autocorrelationMatrixSize = 100 # size of the autocorrelation matrix
    fs=2000 # sampling frequency
    #print(inputDataSamples)
    hankelMatrix=np.zeros((inputDataSamples.size-autocorrelationMatrixSize, 100))
    for i in range(0, inputDataSamples.size-autocorrelationMatrixSize): 
        for j in range(0, autocorrelationMatrixSize):
            hankelMatrix[i,j]=inputDataSamples[i+j+1]
            j=j+1
        i=i+1
    
    Rx=(1/(inputDataSamples.size-autocorrelationMatrixSize))*((hankelMatrix.conj().transpose()).dot(hankelMatrix))
    #print("Rx")
    #print(Rx)

    ###eigenvalue decomposistion
    e_lambda, e_Vector = LA.eig(Rx) #e_lambda= eigen values , e_Vector=eigen vector
    #print("Eigen Values")
    #print(e_lambda)
    
    ###eigen value sorting in decending order and their respective eigen vectors with preserved indeces
    idx = e_lambda.argsort()[::-1]   #sort eigen values in decending order
    #print("idx")
    #print(idx)
    s_lambda = (e_lambda[idx])    #sorted eigen values
    sorted_E_Vector = e_Vector[:,idx]  #sorted eigen vector
    #print("sorted eigen values")
    #print(s_lambda)
    #print("sorted eigen vectors")
    #print(sorted_E_Vector)

    ### Implementing equation
    RD=np.zeros(math.ceil(autocorrelationMatrixSize/2)-1)
    RDI=np.zeros(math.ceil(autocorrelationMatrixSize/2)-1)
    #print("Rd")
    #print(RD)
    #Remember!! range() function in python doesn't return end value so -1 instead of-2
    for i in range(0, math.ceil(autocorrelationMatrixSize/2)-1): 
        a1 = (s_lambda[2*i]+s_lambda[2*i+1])/2
        a2 = (s_lambda[2*i+2]+s_lambda[2*i+3])/2
        a = (a1-a2)/a1
        RDI[i] = i
        RD[i] = a
    
    #print("Rd")
    #print(RD)
    #print("Rdi")
    #print(RDI)

    ###sorting RD in descending order.
    #sorted_RD=np.sort(RD)[::-1]
    #print("sorted_RD")
    #print(sorted_RD)

    ###finding out then largest RD and corresponding RDI
    maxRDI=np.argmax(np.abs(RD),axis=0)
    maxRD=np.amax(np.abs(RD))
    #print("maxRDI")
    #print(maxRDI)
    #print("maxRD")
    #print(maxRD)
    #ploting RD vs RDI
    plt.plot(np.abs(RD), label='RD')
    plt.xlabel('Rank Index (RDI)')
    plt.ylabel('Rank Difference (RD)')
    plt.savefig('static/images/RD_RDI_plot.png')
    plt.close()
    order=maxRDI+1
    print(order)
        
        ###esprit
    Esig=np.zeros((autocorrelationMatrixSize,2*order))

    for i in range(0, autocorrelationMatrixSize):
        for j in range(0, 2*order):
            Esig[i,j] = sorted_E_Vector[i,j]
            j = j+1 
        i= i+1
        #print("Esig")
        #print(Esig.shape)
    ds = 1
    iM=np.eye(autocorrelationMatrixSize-ds)#
        #print(iM)

        ###creating the distance matrix.
    ds = 1;
    ods=np.zeros((1,autocorrelationMatrixSize-ds ))
        #print(ods)

        ###transposing the distance matrix
    od=np.transpose(ods)
        #print(od.shape)
        ###transposing the distance matrix
    s1=np.concatenate((iM,od),1)
    s2=np.concatenate((od,iM),1)

    R1=s1.dot(Esig)
    R2=s2.dot(Esig)

    Psi1=(R1.conj().transpose()).dot(R1)
        #print("Psi1")
        #print(Psi1)
    Psi1_inv = LA.inv(Psi1)
    Psi2=(R1.conj().transpose()).dot(R2)
    Psi = Psi1_inv.dot(Psi2)
        #print("Psi")
        #print(Psi)

    Z,vec = LA.eig(Psi)
        #print("Z")
        #print(Z)
    temp=complex(1,1)
    temp_frequencies1=np.zeros((2*order,1 ))
    for f in range(0, 2*order):
        temp=(cmath.log(Z[f]))
        temp_frequencies1[f] = (temp.imag)/(2* math.pi)
        f=f+1
        ##new code temp_frequencies1=fre; temp_frequencies=freq
    temp_frequencies=fs*temp_frequencies1
    sorted_frequencies=-np.sort(-temp_frequencies, axis=0) #sort in decendng order
    frequencies=np.zeros((int(sorted_frequencies.size/2),1 ))
    for i in range(0, frequencies.size):
        frequencies[i]=sorted_frequencies[i]
    #frequencies_rounded=np.around(frequencies)

    print("frequencies")
    print(frequencies)
    #print("Rounded frequencies")
    #print(frequencies_rounded)

    frequencies_without_fs=-np.sort(-temp_frequencies1, axis=0) #sort frequencies without fs in decendng order

    Rs1=Rx[:,autocorrelationMatrixSize-2*order:autocorrelationMatrixSize]
        #print(Rs1)

    V=np.zeros((autocorrelationMatrixSize,2*order ),dtype=np.complex128)
        #print(frequencies_without_fs.shape)

    for k in range(0, autocorrelationMatrixSize):
        tempx=np.exp(((2*math.pi*(k))*frequencies_without_fs)*1j)
          #print(tempx)
        V[[k],:]=tempx.reshape(1,-1) #convert tempx column vector to row vector by reshapeing and copy the resultant as a K'th row

    VH=(V.conj().transpose()).conj()
    Amplitude=(LA.inv(VH.dot(V))).dot(VH).dot(Rs1)
    namp=(2*np.sqrt(Amplitude[:,0])).conj().transpose()
    namp1=np.abs(namp)
    amplitudes=np.zeros((int(namp1.size/2),1 ))
    for i in range(0, amplitudes.size):
        amplitudes[i]=namp1[i]
    #frequencies_rounded=np.around(frequencies)
    print("Amplitude")
       #print(namp)
    print(amplitudes)
    frequencies = np.array(frequencies).flatten() 
    amplitudes = np.array(amplitudes).flatten()
    freq_amp_pairs = list(zip(frequencies, amplitudes))
    # Save the plot
    plt.bar(frequencies, amplitudes, color='cornflowerblue', width=3.5)
    plt.xlabel('Frequencies (rounded)')
    plt.ylabel('Amplitude')
    plt.title('Frequency vs Amplitude')
    plt.legend(['Amplitude'])
    plt.grid(True)
    plt.xticks(frequencies, [f"{freq:.2f}" for freq in frequencies])
    plt.savefig('static/images/Frequency_Amplitude_plot4.png')
    plt.close()
    
    # Return the calculated frequencies and amplitudes
    #return render_template("index.html", freq_amp_pairs=freq_amp_pairs)

    # plt.show()

def FFT(Voltage,fs):
     
#     inputSignal = pd.read_csv("abc.csv",header=None) #reading uploaded input csv file
#     inputDataSamples=inputSignal.iloc[0,:]  #array containing inputsignal)
#     data=np.transpose(inputDataSamples)


    Fs = 2000  # Sampling frequency in Hz
    N = len(Voltage)  # Number of samples
    T = 1 / Fs  # Sampling interval in seconds

    # Compute the FFT
    fft_result = np.fft.fft(Voltage)

    # Compute the frequency bins
    freqs = np.fft.fftfreq(N, T)

    # Since FFT output is symmetric, take only the positive frequencies
    positive_freqs = freqs[:N // 2]
    positive_fft = fft_result[:N // 2]

    # Compute the amplitude spectrum
   
    frequencies_rounded = np.round(positive_freqs, decimals=2)
    amplitudes = np.abs(positive_fft) * 2 / N
    #freq_amp_pairs = list(zip(frequencies_rounded, amplitudes))
    # Save the plot
    plt.figure(figsize=(12, 6))
    plt.bar(frequencies_rounded, amplitudes, width=fs/N*0.8, color='cornflowerblue')
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Amplitude')
    plt.title('Frequency Spectrum')
    plt.grid(True)
    plt.xlim(0, fs / 2)
    plt.legend(['Amplitude'])
    plt.tight_layout()
    plt.savefig('static/images/Frequency_Spectrum.png')
    plt.close()
    print(frequencies_rounded)
    print(amplitudes)
    # Return the calculated frequencies and amplitudes
    #return render_template("index.html", freq_amp_pairs=freq_amp_pairs)
@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        # Handle file upload
        file = request.files['file']
        frequency = float(request.form['frequency'])
        
        # Save and read the uploaded file
        if file and file.filename.endswith('.csv'):
            filename = secure_filename(file.filename)
            filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
            file.save(filepath)
            inputSignal = pd.read_csv(filepath, header=None)
            inputDataSamples = inputSignal.iloc[0, :]
            Voltage = np.transpose(inputDataSamples)

            x1 = Voltage  # Assuming Voltage is defined somewhere in your Python environment
            x = x1[:1000]
            fs = int(frequency)
            ts = 1 / fs

            # Time vector
            N = 1000
            t = np.arange(0, N * ts, ts)

            # Frequency and calculations for samples per cycle
            f = 50
            Nc = t[-1] * f  # No of cycles in given window
            Npc = round(N/Nc)

            # Downsampling x
            xp1 = []
            y = np.zeros(len(x))

            for i in range(0, len(t), Npc):
                xp1.append(x[i])

            # Rounding xp1
            xp2 = np.round(xp1, 5)
            xp=xp2[1:len(xp2)]

            # Calculate derivative of x
            for i in range(2, len(t)):
                y[i-1] = (x[i]-x[i-1])/ ts


            # Downsampling y
            yp1 = []
            for i in range(0, len(y), Npc):
                yp1.append(y[i])

            # Rounding yp1
            yp2 = np.round(yp1, 5)
            yp=yp2[1:len(yp2)]

            # Plotting
            plt.plot(xp, yp, 'o', color = 'cornflowerblue')
            plt.xlabel('xp')
            plt.ylabel('yp')
            plt.title('Plot of downsampled xp vs. yp')
            plt.savefig('static/images/xp_vs_yp.png')
            plt.show()
            plt.close()

            # Calculate sum for y
            ys = np.sum(np.abs(yp[1:] - yp[0]))

            # Calculate sum for x

            xs = np.sum(np.abs(xp[1:] - xp[0]))
            #plot_type = "Only harmonics are present"

            # Check for interharmonics
            if abs(ys) + abs(xs) > 0.1:
                print('Interharmonics are present')
                MEMO_ESPRIT(Voltage,fs)
                plot_type = "Interharmonics are present"

                xp_yp_url = url_for('static', filename='images/xp_vs_yp.png')
                rd_rdi_url = url_for('static', filename='images/RD_RDI_plot.png')
                plot_url = url_for('static', filename='images/Frequency_Amplitude_plot4.png')
                return render_template('index.html', plot_type=plot_type, plot_url=plot_url, xp_yp_url=xp_yp_url, rd_rdi_url=rd_rdi_url)

            else:
                FFT(Voltage, fs)
                plot_type = "Only harmonics are present"

                xp_yp_url = url_for('static', filename='images/xp_vs_yp.png')
                plot_url = url_for('static', filename='images/Frequency_Spectrum.png')
                return render_template('index.html', plot_type=plot_type, plot_url=plot_url, xp_yp_url=xp_yp_url)


    return render_template("index.html")

if __name__ == "__main__":
    app.run(debug=True)
