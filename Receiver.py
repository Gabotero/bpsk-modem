# Wrote by Gabriel Otero Pérez
# 2025

import numpy as np
import adi
from multiprocessing import Process, Queue
import signal
from collections import deque
import dearpygui.dearpygui as dpg
import time
from scipy.signal import resample_poly

#delete
import time
import random

class Demod(Process):

	def __init__(self, samples_buffer, matched_filter, downsample_factor, sps, processing_block_size, iq_output_buffer, bits_buffer):
		Process.__init__(self)
		self.samples_buffer = samples_buffer
		self.matched_filter = matched_filter
		self.downsample_factor = downsample_factor
		self.sps = sps
		self.block_size = processing_block_size
		self.iq_output_buffer = iq_output_buffer
		self.bits_buffer = bits_buffer

		## Mueller and Muller parameters
		# To be saved between iterations
		self.mu = 0.0						# Guess of the timing phase offset/error
		self.mu_gain = 0.175					# Proportional part K_p of the PLL
		self.omega = self.sps					# Estimate of the number of samples per symbol
		self.omega_gain = 0.25 * self.mu_gain * self.mu_gain	# Gain setting for the integrator constant K_i (see phased locked loops)
		self.omega_relative_limit = 0.005			# Sets the maximum relative deviation from Omega and depends on the clock specifications, i.e., ppm
		self.omega_mid = self.omega
		self.omega_limit = self.omega_mid * self.omega_relative_limit

		self.orphan_samples = 0				# Number of samples pending to merge with the next block
		self.interpolation_factor = 1				# Not used for now

		# Initialized but need to be restarted
		self.ii = 0 # Input index
		self.oo = 2 # Output index
		self.orphans = np.zeros(0, dtype=complex)

		# Auxiliary vectors to save the inputs that we need in each iteration as well as the decisions we make. Shouldn't need > (block_size + orphan_samples(< sps) ) length. 
		self.inputs = np.zeros(self.block_size + self.sps, dtype=complex)		# Grabbed samples
		self.decisions = np.zeros(self.block_size + self.sps, dtype=complex)	# Sliced symbols


		## AGC parameters
		self.agc_gain = 0.005
		self.e_gain = 1.0
		self.min_e_gain = -20.0
		self.max_e_gain = 20.0
		#2
		self.agc_alpha = 0.01
		self.avg_power = 0
		self.reference = 1.0


		## Costas Loop parameters
		#To be saved between iterations
		self.costas_alpha = 0.132
		self.costas_beta = 0.00932
		self.costas_phase_offset = 0.0
		self.costas_freq_offset = 0.0

	def branchless_clip(self, x, clip):
		return 0.5 * (abs(x + clip) - abs(x - clip))

	def RMS_AGC(self, samples):

		output = np.zeros_like(samples)
		for i in range(len(samples)):
			sample = samples[i]*self.e_gain
			#sample = np.clip(sample, -1e10, 1e10)
			self.e_gain += (1-np.abs(sample))*self.agc_gain
			#self.e_gain = max(self.min_e_gain, min(self.e_gain, self.max_e_gain))
			output[i] = sample

		return output

	def RMS_AGC2(self, samples):

		output = np.zeros_like(samples)
		for i in range(len(samples)):

			self.avg_power = (1-self.agc_alpha)*self.avg_power + self.agc_alpha*(np.abs(samples[i])**2)
			rms = np.sqrt(self.avg_power)
			output[i] = (samples[i]/rms)*self.reference

		return output

	def MM_Sync(self, samples):

		# Interpolate?
		#samples = signal.resample_poly(samples, UP, DOWN) # Check here for GNU Radio implementation of interpolation M&M https://edfuentetaja.github.io/sdr/m_m_gnu_radio_analysis/. Also, note the class variable self.interpolation_factor to adjust it, if needed.

		samples = np.concatenate((self.orphans, samples))

		#print(f"Processing {len(samples)} samples...")

		# Reset the auxiliary vectors
		self.inputs = np.concatenate((self.inputs[self.oo-2:self.oo], np.zeros(len(samples), dtype=complex)))
		self.decisions = np.concatenate((self.decisions[self.oo-2:self.oo], np.zeros(len(samples), dtype=complex)))

		# Reset the output index. First 2 samples come from the previous block
		self.oo = 2

		num_inputs = len(samples) - self.sps # If TX transmits 9600 symbols/s and we sample at (1e6/downsample) = (1e6/8) = (125e3), then we have (125e3/9600) ~ 13 samples per symbol (self.sps). Therefore, process samples until I don't have enough for another symbol, i.e., to "jump" to the next symbol.

		while(self.ii < num_inputs and self.oo < len(self.inputs)):
			# New input to the synchronizer
			self.inputs[self.oo] = samples[self.ii + int(self.mu)] # If interpolating: samples[self.ii*self.interpolation_factor + int(self.mu*self.interpolation_factor)]

			# Slicer
			self.decisions[self.oo] = int(np.real(self.inputs[self.oo]) > 0) + 1j * int(np.imag(self.inputs[self.oo]) > 0)

			# Compute the error
			y = (self.inputs[self.oo] - self.inputs[self.oo-2]) * np.conj(self.decisions[self.oo-1])
			x = (self.decisions[self.oo] - self.decisions[self.oo-2]) * np.conj(self.inputs[self.oo-1])
			timing_error = np.real(y-x)

			# Update M&M variables
			self.omega = self.omega + self.omega_gain * timing_error
			self.omega = self.omega_mid + self.branchless_clip(self.omega-self.omega_mid, self.omega_limit)
			self.mu = self.mu + self.omega + self.mu_gain * timing_error
			#print(f"Omega is {self.omega}")

			self.ii += int(np.floor(self.mu))
			self.mu = self.mu - np.floor(self.mu)

			# Saw this clamp in satdump's implementation
			if (self.ii < 0):
				self.ii = 0

			self.oo += 1

			#print(f"mu is: {self.mu}. Omega is: {self.omega}")

		#print(f"Termino y mu es {self.mu}")
		# Check whether we need orphan samples for the next block
		if(self.ii >= len(samples)): # We don't need the remaining samples and we need to jump some of the new
			self.ii = self.ii - len(samples)
			self.orphans = np.zeros(0, dtype=complex)

		else:	# We do need the orphan samples so that we don't lose continuity
			self.orphans = samples[self.ii:]
			self.ii = 0


		#print(f"mu is: {self.mu}. Omega is: {self.omega}")
		# Return the output
		return self.inputs[2:self.oo] # self.oo is updated but not included


	def Costas_Loop(self, samples):

		num_samples = len(samples)
		output = np.zeros(num_samples, dtype=complex)

		for i in range(num_samples):

			# Correct the input sample according to the estimated phase offset
			output[i] = samples[i] * np.exp(-1j*self.costas_phase_offset)

			# Compute the new error estimate
			error = self.error_detector(output[i])

			# Update the loop
			self.costas_freq_offset += self.costas_beta * error
			self.costas_phase_offset += self.costas_freq_offset + (self.costas_alpha * error)

			while self.costas_phase_offset >= 2*np.pi:
				self.costas_phase_offset -= 2*np.pi
			while self.costas_phase_offset < 0:
				self.costas_phase_offset += 2*np.pi

		return output


		#TODO: output frequency offset evolution to a buffer for plotting.


	def error_detector(self, sample): #BPSK

		err = np.real(sample) * np.imag(sample)
		return err

	def slicer(self, complex_samples): # Only BPSK for now.

		# If BPSK...
		real_parts = np.real(complex_samples)
		bits = np.where(real_parts >= 0, 1, 0)
		# else...

		return bits


	def run(self):

		while(1):

			iq_data = self.samples_buffer.get() # Gets a block of size "FFT size"

			# Coarse Frequency Sync (i.e., Doppler) TODO

			# Decimate
			iq_data_dec = resample_poly(iq_data, 1, self.downsample_factor) # signal.resample_poly(x, up, down, window, padtype)d
			# Sample rate here should be 1 MHz/8 = 125 KHz and we should have ~13 samples per symbol (self.sps)

			# Apply matched filter
			#iq_data_filtered = np.convolve(iq_data_dec, self.matched_filter, 'same') #TODO: Check {‘full’, ‘valid’, ‘same’} if needed. Default is Full

			# AGC
			data_controlled_gain = self.RMS_AGC2(iq_data_dec)

			# Clock Sync (Mueller and Muller)
			data_time_synced = self.MM_Sync(data_controlled_gain)

			# Fine Frequency Sync (Costas Loop)
			data_freq_synced = self.Costas_Loop(data_time_synced)

			# Output of demod
			self.iq_output_buffer.put(data_freq_synced)
			self.bits_buffer.put(self.slicer(data_freq_synced))

class BitStreamProcessing(Process):

	def __init__(self, bitstream_buffer, syncWord):
		Process.__init__(self)

		self.bitstream = bitstream_buffer
		self.syncWord = int(syncWord, 16)
		self.syncLength = 64
		self.syncWordInverted = ~self.syncWord & ((1 << self.syncLength) - 1) # 0xFCB88938D8D76A4F
		self.buffer_size = self.syncLength
		self.data_register = 0

		self.state = "SEARCHING_SYNC"
		self.flip_bits = 0 # To remove abiguity
		self.bit_counter = 0

		self.frame_size = 96 # [bits]
		self.frame_buffer = []

		self.state_switch = {
			'SEARCHING_SYNC' : self.searching_sync,
			'GATHERING_TM_TRANSFER_FRAME' : self.gathering_TM_Transfer_frame
		}



	def searching_sync(self):

		if(self.data_register == self.syncWord):
			print("Sync detected!")
			self.state = 'GATHERING_TM_TRANSFER_FRAME'
		elif(self.data_register == self.syncWordInverted):
			print("Inverted Sync detected! Flipping bits...")
			self.state = 'GATHERING_TM_TRANSFER_FRAME'
			self.flip_bits ^=1 # Invert bits if an inverted sync was received.

	def gathering_TM_Transfer_frame(self):

		self.bit_counter += 1

		self.frame_buffer.append(self.data_register & 0x1)

		if(self.bit_counter == self.frame_size): # TODO: send raw frame to another process/machine for further decoding.

			chars = []
			for b in range(0, len(self.frame_buffer), 8):  # Process bytes
				byte = self.frame_buffer[b:b+8]
				byte_str = ''.join(str(bit) for bit in byte)
				chars.append(chr(int(byte_str, 2)))
			print("The received text is: ", ''.join(chars))

			self.state = 'SEARCHING_SYNC'
			self.bit_counter = 0
			self.frame_buffer = []

	def run(self):

		while(1):

			while(self.bitstream.qsize() > 0):

				new_bits = self.bitstream.get()

				for bit in new_bits:

					bit = int(bit) ^ self.flip_bits #  flip_bits will tell whether we flip the bits or not for reception.
					# DBPSK
					#received_bit = int(bit) ^ self.last_bit
					#self.last_bit = received_bit


					# Store a new bit
					self.data_register = ((self.data_register << 1) | (int(bit) & 0x1)) & ((1 << self.buffer_size) - 1)

					# Process the new data
					self.state_switch[self.state]()



class PlutoSDR(Process):

	def __init__(self, sample_rate, center_freq, num_samps, hw_gain, buff, buffer_fft):
		Process.__init__(self)
		self.sample_rate = sample_rate
		self.center_freq = center_freq
		self.num_samps = num_samps # number of samples we need per FFT
		self.hw_gain = hw_gain
		self.buff = buff
		self.buffer_fft = buffer_fft

	def setup(self):
		sdr = adi.Pluto('ip:192.168.2.1')
		sdr.gain_control_mode_chan0 = 'manual'
		sdr.rx_hardwaregain_chan0 = self.hw_gain # dB
		sdr.rx_lo = int(self.center_freq)
		sdr.sample_rate = int(self.sample_rate)
		sdr.rx_rf_bandwidth = int(self.sample_rate) # filter width, just set it to the same as sample rate for now
		sdr.rx_buffer_size = 32768 #self.num_samps # number of samples returned per call to rx()

		return sdr

	def run(self):

		sdr = self.setup()
		buffer_samples = np.array([], dtype=np.complex64)

		while(1):

			# Receive a new block of samples (32768)
			samples = sdr.rx()

			# Fill in the buffer
			buffer_samples = np.concatenate((buffer_samples, samples))


			while len(buffer_samples) >= self.num_samps:
				block = buffer_samples[:self.num_samps]
				self.buff.put(block)
				self.buffer_fft.put(block)
				buffer_samples = buffer_samples[self.num_samps:]
				#print("Buffer has ", self.buff.qsize(), "blocks")
				#print("Buffer fft has ", self.buffer_fft.qsize(), "blocks")

class FFT_Process(Process):

	def __init__(self, fft_size, num_avg, sample_rate, buff, output_buffer):
		Process.__init__(self)
		self.fft_size = fft_size
		self.sample_rate = sample_rate
		self.numAvg = num_avg
		self.buff = buff
		self.output_buffer = output_buffer

	def run(self):

		PSD_shifted = np.zeros(self.fft_size, dtype=np.float64)
		fft_count = 0
		while(1):
			samples = self.buff.get() # Will block if buffer is empty. Gets fft_size samples.
			samples = samples*np.hamming(self.fft_size) # Apply a Hamming window
			PSD = np.abs(np.fft.fft(samples))**2 / (self.fft_size*self.sample_rate)
			PSD_log = 10.0*np.log10(PSD)
			PSD_shifted += np.fft.fftshift(PSD_log)
			fft_count += 1

			if(fft_count == self.numAvg):
				self.output_buffer.put(PSD_shifted/self.numAvg)
				PSD_shifted = np.zeros_like(PSD_shifted)
				fft_count = 0

			#print("We have ", self.output_buffer.qsize(), " FFTs")



def generate_matched_filter(samps_per_symbol, n_taps = 101, beta = 0.35):

	Ts = samps_per_symbol # Assume sample rate is 1 Hz, so sample period is 1s, so the symbol period is equal to the samples per symbol
	t = np.arange(n_taps) - (n_taps-1)//2
	h = (1/Ts)*np.sinc(t/Ts) * np.cos(np.pi*beta*t/Ts) / (1 - (2*beta*t/Ts)**2)

	return h


if __name__ == "__main__":

	samp_rate = 1e6
	center_frequency = 145e6
	fftSize = 4096
	num_avg_ffts = 20.0 # For visualization purposes.
	hardwareGain = 60.0


	# BPSK parameters
	baudrate = 9600
	target_samples_per_symbol = 13
	downsample = int(samp_rate/(baudrate*target_samples_per_symbol))
	print("Downsample is now: ", str(downsample))
	mtchd_filter = generate_matched_filter(samps_per_symbol = target_samples_per_symbol, n_taps = 101, beta = 0.35) # We can assume a sample period of 1 second to “normalize”. It means our symbol period Ts is 13 because we have 13 samples per symbol
	syncWord_hex = "034776C7272895B0"

	sampling_buffer = Queue()
	sampling_buffer_fft = Queue()
	fft_output_buffer = Queue()
	mm_iq_output_buffer = Queue()
	binary_buffer = Queue()


	samplerProcess = PlutoSDR(sample_rate = samp_rate, center_freq = center_frequency, num_samps = fftSize, hw_gain = hardwareGain, buff = sampling_buffer, buffer_fft = sampling_buffer_fft)
	fftProcess = FFT_Process(fft_size = fftSize, num_avg = num_avg_ffts, sample_rate = samp_rate, buff = sampling_buffer_fft, output_buffer = fft_output_buffer)
	demodProcess = Demod(samples_buffer = sampling_buffer, matched_filter = mtchd_filter, downsample_factor = downsample, sps = target_samples_per_symbol, processing_block_size = fftSize, iq_output_buffer = mm_iq_output_buffer, bits_buffer = binary_buffer)
	bitStreamProcess = BitStreamProcessing(bitstream_buffer = binary_buffer, syncWord = syncWord_hex)

	samplerProcess.start()
	fftProcess.start()
	demodProcess.start()
	bitStreamProcess.start()

	fft_xAxis = np.fft.fftshift(np.fft.fftfreq(fftSize, 1/samp_rate) + center_frequency)/1e6 # (size, timestep)

	constellation_points_display_limit = 2048
	constellation_plot_queue = deque(maxlen=constellation_points_display_limit)

	dpg.create_context()
	dpg.create_viewport(title='BPSK Modem - 2025', width=1440, height=950)
	dpg.setup_dearpygui()


	with dpg.window(label = "Receiver", width = 1600, height = 1200):
		dpg.add_text("Sample rate is: %.2f [MS/s]. Averaging %d FFTs"%(samp_rate/1e6, num_avg_ffts), tag = "performance_info")
		dpg.add_text("Fps: 000.00f", tag="fps_tag")

		with dpg.theme() as scatter_theme:
			with dpg.theme_component(dpg.mvScatterSeries):
				dpg.add_theme_color(dpg.mvPlotCol_Line, (255, 0, 0), category=dpg.mvThemeCat_Plots)
				#dpg.add_theme_style(dpg.mvPlotStyleVar_Marker, dpg.mvPlotMarker_Square, category=dpg.mvThemeCat_Plots)
				dpg.add_theme_style(dpg.mvPlotStyleVar_MarkerSize, 2, category=dpg.mvThemeCat_Plots)

		with dpg.plot(label = "Spectrum", width = 1024, height = 350, pos=(100, 100)):
			dpg.add_plot_legend()
			dpg.add_plot_axis(dpg.mvXAxis, label = "Frequency [MHz]")
			#dpg.set_axis_limits(dpg.last_item(), (center_frequency - samp_rate/2)/1e6, (center_frequency + samp_rate/2)/1e6)
			dpg.set_axis_limits(dpg.last_item(), (center_frequency - samp_rate/16)/1e6, (center_frequency + samp_rate/16)/1e6)

			dpg.add_plot_axis(dpg.mvYAxis, label = "PSD", tag = "fft_plot_y_axis")
			dpg.set_axis_limits(dpg.last_item(), -50, 50)

			dpg.add_line_series([1, 2, 3, 4], [1, 2, 3, 4], parent="fft_plot_y_axis", tag = "fft_plot")

		with dpg.plot(label = "Constellation diagram", width = 400, height = 400, pos=(0, 100+400)):
			dpg.add_plot_axis(dpg.mvXAxis, label = "In-Phase (I)")
			dpg.set_axis_limits(dpg.last_item(), -2, 2)
			dpg.add_plot_axis(dpg.mvYAxis, label = "Quadrature (Q)", tag = "iq_constellation_y_axis")
			dpg.set_axis_limits(dpg.last_item(), -2, 2)

			scatter = dpg.add_scatter_series([-1, -1, 1, 1], [-1, 1, -1, 1], label="IQ data", parent="iq_constellation_y_axis", tag="iq_plot")

		with dpg.plot(label = "Eye diagram", width = 800, height = 400, pos=(400, 100+400)):
			xaxis = dpg.add_plot_axis(dpg.mvXAxis, label = "Time [samples]")
			dpg.set_axis_limits(dpg.last_item(), 0, 31)
			dpg.add_plot_axis(dpg.mvYAxis, label="Amplitude", tag="eye_diagram_y_axis")
			dpg.set_axis_limits(dpg.last_item(), -2, 2)
			dpg.add_line_series([], [], parent="eye_diagram_y_axis", tag="eye_plot")

	dpg.bind_item_theme(scatter, scatter_theme)

	dpg.show_viewport()

	rendering_time = 0
	frame_counter = 0
	start_time = time.time()

	while dpg.is_dearpygui_running():

		# This runs in the render loop
		# Can be manually stopped by using stop_dearpygui()

		#############
		# Code here #
		#############

		# Update FFT
		if(not fft_output_buffer.empty()): # TODO: test without this double check !=0 and >0
			#print("We have ", fft_output_buffer.qsize(), " FFTs left.")
			while(fft_output_buffer.qsize() > 0):
				dpg.set_value("fft_plot", [fft_xAxis, fft_output_buffer.get()])

		# Update IQ and eye plot
		if(not mm_iq_output_buffer.empty()): # TODO: test without this double check !=0 and >0

			while(mm_iq_output_buffer.qsize() > 0):
				constellation_plot_queue.extend(mm_iq_output_buffer.get())

			#print("We have ", len(constellation_plot_queue), " samples.")

			x_data = [np.real(sample) for sample in constellation_plot_queue]
			y_data = [np.imag(sample) for sample in constellation_plot_queue]

			dpg.set_value("iq_plot", [x_data, y_data])

			#Upsample the real part for the eye plot
			upsampled = resample_poly(x_data, 16, 1)
			dpg.set_value("eye_plot", [np.arange(len(upsampled))%32, np.ascontiguousarray(upsampled)])


		#############


		dpg.render_dearpygui_frame()

		elapsed = time.time() - start_time
		rendering_time += elapsed
		frame_counter+=1

		if(frame_counter == 100):
			fps = 1/(rendering_time/frame_counter)
			dpg.set_value("fps_tag", "Fps: %.2f"%fps)
			rendering_time = 0
			frame_counter = 0

		start_time = time.time()


	dpg.destroy_context()

	bitStreamProcess.terminate()
	demodProcess.terminate()
	fftProcess.terminate()
	samplerProcess.terminate()

