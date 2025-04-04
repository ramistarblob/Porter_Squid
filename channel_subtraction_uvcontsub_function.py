import numpy as np
import argparse

def channel_continuum_subtract(ms_set, spw, ref_freq_GHz):
    #Once spw known, create casa table (ms set name will change) CASA TABLE DOCUMENTATION
    tb.open(ms_set +'/SPECTRAL_WINDOW')

    #REF FREQ OF OUR ISOTOPOLOQUE (will change) [in hz]
    iso_ref_freq =  ref_freq_GHz*10**9

    #channel frequencies, row # = SPW (will change)
    chan_freqs = tb.getcol(columnname="CHAN_FREQ", startrow=spw, nrow=1, rowincr=1)
	
    #central freq
    cent_freq = chan_freqs[4095]
    print('CEN FREQ =', cent_freq*10**-9)
	
    #CHAN WIDTH row # = SPW (will change)
    chan_width = tb.getcol(columnname="CHAN_WIDTH", startrow=spw, nrow=1, rowincr=1)[0]
    print('CHAN Width (GHz) =', chan_width*10**-9)
    
    #CHANNEL 0 frequency
    #chan0_freq = chan_freqs[0]
    #print('CHAN 0 FREQ (GHz) =', chan0_freq*10**-9)
    chan0_freq = cent_freq-4095*chan_width
    print('CHAN 0 FREQ (GHz)', chan0_freq*10**-9)
	
    chan_last_freq = chan0_freq + (8191 * chan_width)
    print('Last Channel Frequency (GHz) [8191] =', chan_last_freq * 10**-9)


    iso_corres_channel = (((iso_ref_freq-chan0_freq)/(chan_width)))
    print('Isotopologue corres channel',np.round(iso_corres_channel,decimals=0))
    #print('chans to exc',np.round(iso_corres_channel,decimals=0)-51)
    #print('chans to exc freq', (np.round(iso_corres_channel,decimals=0)-51)*10**-9+chan0_freq*10**-9)
    #print('chans to exc',np.round(iso_corres_channel,decimals=0)+51)
    #print('chans to exc freq', (np.round(iso_corres_channel,decimals=0)+51)*10**-9+chan0_freq*10**-9)

    line_100mhz_chann_a = np.round((iso_ref_freq-chan0_freq-1*10**8)/(chan_width))
    print(line_100mhz_chann_a)
    #print(chan0_freq*10**-9+line_100mhz_chann_a*chan_width*10**-9)
    
    line_100mhz_chann_b = np.round((iso_ref_freq-chan0_freq+1*10**8)/(chan_width))
    print(line_100mhz_chann_b)
    #print(chan0_freq*10**-9+line_100mhz_chann_b*chan_width*10**-9)
	
    print(f"Emission may be present from chan {line_100mhz_chann_a} ({chan0_freq*10**-9+line_100mhz_chann_a*chan_width*10**-9} GHz) to {line_100mhz_chann_b} ({chan0_freq*10**-9+line_100mhz_chann_b*chan_width*10**-9} GHz)")
    
    if line_100mhz_chann_a < line_100mhz_chann_b:

    	print('Channels to subtract', f'20~{(line_100mhz_chann_a)-1};{(line_100mhz_chann_b)+1}~{len(chan_freqs)-20}')
    else:
    	print('Channels to subtract', f'20~{(line_100mhz_chann_b)-1};{(line_100mhz_chann_a)+1}~{len(chan_freqs)-20}')


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('ms_set')
	parser.add_argument('spw', type=int)
	parser.add_argument('ref_freq_GHz', type=float)
	args = parser.parse_args()
	channel_continuum_subtract(args.ms_set, args.spw, args.ref_freq_GHz)



