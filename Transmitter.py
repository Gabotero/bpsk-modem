import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy.signal import resample_poly
from matplotlib.animation import FuncAnimation
import adi
import time

def generate_matched_filter_tx(samps_per_symbol, n_taps = 101, beta = 0.35):

	Ts = samps_per_symbol # Assume sample rate is 1 Hz, so sample period is 1s, so the symbol period is equal to the samples per symbol
	t = np.arange(n_taps) - (n_taps-1)//2
	h = (1/Ts)*np.sinc(t/Ts) * np.cos(np.pi*beta*t/Ts) / (1 - (2*beta*t/Ts)**2)
	h*=Ts
	return h

# Text to be transmitted
message = "Hi there! :)"
message_bits = ''.join(format(ord(c), '08b') for c in message)
message_bits = [int(bit) for bit in message_bits]  # Convertir a lista de enteros

print("Message size is: ", len(message_bits))

# Sync word
syncWord_hex = "034776C7272895B0"
syncWord_bits = ''.join(format(int(syncWord_hex[i:i+2], 16), '08b') for i in range(0, len(syncWord_hex), 2))
syncWord_bits = [int(bit) for bit in syncWord_bits]

# Combine both
bits = syncWord_bits + message_bits
num_symbols = len(bits)

print("num symbols is", num_symbols)
samps_per_symbol = 8
baudrate = 9600

signal = np.array([])

for bit in bits:
	pulse = np.zeros(samps_per_symbol) #np.ones(samps_per_symbol) * (bit * 2 -1)
	pulse[0] = bit * 2 - 1 # set the first value to either a 1 or -1
	signal = np.concatenate((signal, pulse)) # add the 8 samples to the signal

#plt.figure(0)
#plt.plot(signal, '.-')
#plt.grid(True)
#plt.show()

n_taps = 101
beta = 0.35
shape_filter_taps = generate_matched_filter_tx(samps_per_symbol, n_taps, beta)

#plt.figure(1)
#plt.plot(np.arange(n_taps) - (n_taps-1)//2, shape_filter_taps, '.')
#plt.grid(True)
#plt.show()


signal_shaped = np.convolve(signal, shape_filter_taps, 'same')

#plt.figure(2)
#plt.plot(signal_shaped, '.-')
#for i in range(num_symbols):
#    plt.plot([i*samps_per_symbol+n_taps//2,i*samps_per_symbol+n_taps//2], [0, signal_shaped[i*samps_per_symbol+n_taps//2]])
#plt.grid(True)
#plt.show()


signal_shaped_iq = signal_shaped + 0.2j * np.ones(len(signal_shaped))

sdr_samp_rate = 1e6

upsample_factor = int(sdr_samp_rate/(samps_per_symbol*baudrate))

print("Upsampling by: ", upsample_factor)



upsampled = resample_poly(signal_shaped_iq, upsample_factor, 1)
upsampled *= 2**14

center_freq = 145e6 # Hz
sdr = adi.Pluto("ip:192.168.4.1")
sdr.sample_rate = int(sdr_samp_rate)
sdr.tx_rf_bandwidth = int(sdr_samp_rate) # filter cutoff, just set it to the same as sample rate
sdr.tx_lo = int(center_freq)
sdr.tx_hardwaregain_chan0 = 0 # Increase to increase tx power, valid range is -90 to 0 dB


sdr.tx_cyclic_buffer = True

# Send data cyclically
print("Transmitting...")
sdr.tx(upsampled)
print("Finished")
time.sleep(10000)

print("Transmission Completed, closing")
