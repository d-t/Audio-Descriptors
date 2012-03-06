/**
 * ##library.name##
 * ##library.sentence##
 * ##library.url##
 *
 * Copyright ##copyright## ##author##
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General
 * Public License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place, Suite 330,
 * Boston, MA  02111-1307  USA
 * 
 * @author      ##author##
 * @modified    ##date##
 * @version     ##library.prettyVersion## (##library.version##)
 */

package audio.descriptors;


import java.util.ArrayList;

import ddf.minim.AudioSource;
import ddf.minim.analysis.FFT;
import processing.core.*;

/**
 * <code>AudioDescriptor</code> is the class that contains 
 * all the functions that provide low-level audio signal features. 
 * These features give useful information when you want to analyze 
 * a digital audio signal.
 * 
 *  Audio Descriptors can be used for many purposes, such as:
 *  - Speech Recognition
 *  - Music Information Retrieval applications
 *  - Mood Detection
 *  Ð Event Detection
 *  - Pitch Tracking
 *  
 *  This class works with Minim <code>AudioSource</code> objects, 
 *  i.e. <code>AudioInput</code>, <code>AudioOutput</code>, <code>AudioPlayer</code> 
 *  and <code>AudioSample</code>.
 *
 */

public class AudioDescriptor {
	
	// myParent is a reference to the parent sketch
	PApplet myParent;
	
	public final static String VERSION = "##library.prettyVersion##";
	

	/**
	 * Default constructor
	 * 
	 * @example Hello
	 * @param theParent
	 */
	public AudioDescriptor(PApplet theParent) {
		myParent = theParent;
		welcome();
	}
	
	/**
	 * No-parameter constructor
	 */
	public AudioDescriptor(){
		
	}
	
	
	private void welcome() {
		System.out.println("##library.name## ##library.prettyVersion## by ##author##");
	}
	
	
	
	/**
	 * return the version of the library.
	 * 
	 * @return String
	 */
	public static String version() {
		return VERSION;
	}
	
	
	/**
	 * The spectral brightness is defined as the 
	 * amount of energy in the spectrum starting 
	 * from a specified cutoff frequency 
	 * (in this case 1500 Hz), divided by the 
	 * total energy of the signal spectrum
	 * 
	 *  @param x
	 * 			The time-domain signal to be analyzed
	 *  @return b
	 *  		Spectral brightness coefficient
	 */

	public float getRightBrightness(AudioSource x){
	  // compute FFT
	  FFT f = new FFT(x.bufferSize(), x.sampleRate());
	  f.window(FFT.HAMMING);
	  f.forward(x.right);
	  int ss = f.specSize();  // spectrum size
	  
	  final int fc = 1500;  // Why 1500 Hz? By definition...
	  final int kc = f.freqToIndex(fc);
	  
	  float num = 0;  // numerator
	  float den = 0;  // denominator
	  float b = 0;    // result
	  
	  for(int i=0; i<ss; i++){
	    den += f.getBand(i);
	    if(i >= kc){
	      num += f.getBand(i);
	    }
	  }
	  
	  b = num / den;
	  return(b);
	}//end getRightBrightness
	
	
	/**
	 * The spectral centroid is the baricenter 
	 * of the spectrum; it is calculated as the 
	 * weighted mean of the frequencies present 
	 * in the signal, with their magnitude 
	 * as the weights.
	 * 
	 * The result is expressed in Hertz [Hz]
	 * 
	 *   @param x
	 * 			The time-domain signal to be analyzed
	 *   @return asc
	 * 			Audio Spectral Centroid measured in Hertz [Hz]
	 */

	public float getRightCentroid(AudioSource x){
	  // compute FFT
	  FFT f = new FFT(x.bufferSize(), x.sampleRate());
	  f.window(FFT.HAMMING);
	  f.forward(x.right);
	  int Nft = f.specSize();  // spectrum size
	  
	  float deltaF = x.sampleRate() / Nft;  // frequency interval between two FFT bins
	  int kLow = (int) Math.floor(62.5 / deltaF);  // preventing disproportionate weight below this frequency
	  
	  float asc = 0;  // result
	  float num = 0;
	  float den = 0;
	  
	  // compute the first component
	  float Pk = 0;
	  if(kLow > 0){
	    for(int i=0; i<kLow; i++){  // for each band
	      Pk += f.getBand(i);
	    }
	    num += 31.25 * Pk;  // 31.25 Hz = fLow
	    den += Pk;
	  }
	  
	  // compute the remaining components
	  float fk = 0;
	  for(int i=kLow; i<Nft/2; i++){
	    Pk = f.getBand(i);
	    fk = f.indexToFreq(i);
	    num += fk * Pk;
	    den += Pk;
	  }
	  
	  asc = num / den;
	  return(asc);
	  
	}//end getRightCentroid
	
	
	/**
	 * The Spectral Entropy is computed as the 
	 * ShannonÕs entropy of the signal spectrum.
	 * 
	 * The entropy is divided by the logarithm of 
	 * the length N of the spectrum in order to have 
	 * a result which is independent from the windowing length.
	 * 
	 *   @param x
	 * 			The time-domain signal to be analyzed
	 * 	 @return ent
	 * 			Spectral Entropy coefficient
	 */

	public float getRightEntropy(AudioSource x){
	  // compute FFT
	  FFT f = new FFT(x.bufferSize(), x.sampleRate());
	  f.window(FFT.HAMMING);
	  f.forward(x.right);
	  int ss = f.specSize();  // spectrum size
	  
	  float ent = 0;  // result
	  float num = 0;
	  float den = 0;
	  
	  float Xi = 0;
	  for(int i=0; i<ss; i++){
	    Xi = f.getBand(i);
	    num += (Xi * Math.log(Xi));
	  }
	  
	  den = (float) Math.log(ss);  // logarithm of the spectrum length
	  
	  ent = - (num / den);
	  return(ent);
	}//end getRightEntropy
	
	
	/**
	 * Spectral flatness provides a way to quantify 
	 * how tone-like a sound is, as opposed to 
	 * being noise-like. The meaning of tonal in this 
	 * context is in the sense of the amount of peaks 
	 * or resonant structure in a power spectrum, as 
	 * opposed to flat spectrum of a white noise.
	 * 
	 *   @param x
	 * 			The time-domain signal to be analyzed
	 *	 @return asf
	 * 			Spectral Flatness coefficient
	 */

	public float getRightFlatness(AudioSource x){
	  // compute FFT
	  FFT f = new FFT(x.bufferSize(), x.sampleRate());
	  f.window(FFT.HAMMING);
	  f.forward(x.right);
	  
	  float num = 1;
	  float den = 0;
	  float Si = 0;
	  float asf = 0;  // result

	  final int n = -8;
	  final int B = 24;  // number of bands
	  final float loF = (float) (Math.pow(2, n/4.0) * 1000);  // lowest frequency [Hz]
	  final float hiF = (float) (Math.pow(2, B/4.0) * loF);  // highest frequency [Hz]
	  final int loK = f.freqToIndex(loF);
	  final int hiK = f.freqToIndex(hiF);
	  final float reduceFactor = hiK - loK + 1;
	  
	  for(int k=loK; k<hiK; k++){
	    Si = f.getBand(k);
	    num *= Math.pow(Si, 1.0/reduceFactor);  // n-th root
	    den += Si;
	  }
	  
	  den /= reduceFactor;
	  asf = num / den;
	  return(asf);
	}//end getRightFlatness
	
	
	/**
	 * Spectral flux is a measure of how quickly 
	 * the power spectrum of a signal is changing, 
	 * calculated by comparing the power spectrum 
	 * of half a frame against the power spectrum 
	 * of the other half.
	 * 
	 *   @param x
	 * 			The time-domain signal to be analyzed
	 * @return flux
	 * 			Spectral Flux coefficient
	 */
	public float getRightFlux(AudioSource x){
	  final int bufSize = x.bufferSize();
	  float[] ab = new float[bufSize];
	  float[] ab1, ab2;
	  ab1 = new float[bufSize / 2];
	  ab2 = new float[bufSize / 2];
	  
	  // split the audio buffer into 2 halves
	  ab = x.right.toArray();
	  for(int i=0; i<(bufSize / 2); i++){
	    ab1[i] = ab[i];
	    ab2[i] = ab[(bufSize / 2) + i];
	  }
	  
	  // compute the 2 spectra
	  FFT f1 = new FFT(ab1.length, x.sampleRate());
	  FFT f2 = new FFT(ab2.length, x.sampleRate());
	  f1.window(FFT.HAMMING);
	  f1.forward(ab1);
	  f2.window(FFT.HAMMING);
	  f2. forward(ab2);
	  int ss = f1.specSize();  // spectrum size
	  
	  // compute the spectral flux
	  float S1 = 0;
	  float S2 = 0;
	  for(int i=0; i<ss; i++){
	    S1 += Math.pow(f1.getBand(i), 2);
	    S2 += Math.pow(f2.getBand(i), 2);
	  }
	  
	  float flux = (float) Math.sqrt(Math.pow(S2 - S1, 2));
	  return(flux);
	}//end getRightFlux
	
	
	/**
	 * The harmonic ratio (HR) is a measure of the proportion
	 * of harmonic components in the power spectrum. It is in
	 * practice a coefficient associated to each signal frame,
	 * computed as the maximum value of the autocorrelation 
	 * function.
	 * 
	 *   @param x
	 * 			The time-domain signal to be analyzed
	 * 	 @return hr
	 * 			Harmonic Ratio coefficient
	 */

	public float getRightHR(AudioSource x){
	  final int L = x.bufferSize();
	  final int Mmin = 5;  // prevent giving higher importance to low lags
	  final int Mmax = L -1;  // maximum lag
	  
	  float hr = -1;  // result
	  float gamma;  // current autocorrelation value
	  // compute autocorrelation
	  for(int m=Mmin; m<=Mmax; m++){
	    gamma = rightAutoCorr(m, x);
	    if(gamma > hr){
	      hr = gamma;
	    }
	  }
	  
	  return(hr);
	}//end getRightZCR

	/**
	 * Autocorrelation function
	 * 
	 *   @param m
	 *   		Imput lag
	 *   @param x
	 * 			The time-domain signal to be analyzed
	 * @return gamma
	 * 			<code>m</code>-lag autocorrelation output
	 * 
	 * @see <a href="http://en.wikipedia.org/wiki/Autocorrelation">Definition of autocorrelation</a>
	 */
	private float rightAutoCorr(int m, AudioSource x){
	  final int L = x.bufferSize();
	  float num = 0;
	  float den1 = 0;
	  float den2 = 0;
	  float gamma = 0;  // result
	  
	  // compute autocorrelation
	  for(int t=0; t<L-1; t++){
	    if(t+m < L){
	      num += x.right.get(t) * x.right.get(t+m);
	      den1 += Math.pow(x.right.get(t), 2);
	      den2 += Math.pow(x.right.get(t+m), 2);
	    }
	  }
	  
	  gamma = (float) (num / Math.sqrt(den1 * den2));
	  return(gamma);
	}//end autoCorr
	
	
	/**
	 * Inharmonicity estimates the amount of partials that 
	 * are not multiples of the fundamental frequency. It is 
	 * computed by considering the ideal locations of 
	 * harmonics vs the actual harmonics in a spectrum.
	 * 
	 * In this implementation, a lower value implies a more
	 * harmonic signal.
	 * 
	 *   @param x
	 * 			The time-domain signal to be analyzed
	 *	 @return inh
	 *			Inharmonicity coefficient
	 */
	public float getRightInharmonicity(AudioSource x){
	  // compute FFT
	  FFT f = new FFT(x.bufferSize(), x.sampleRate());
	  f.window(FFT.HAMMING);
	  f.forward(x.right);
	  
	  // compute the fundamental frequency f0
	  float f0 = getRightFundFreq(x);
	  
	  // compute the spectrum peaks
	  ArrayList<Integer> peaks = new ArrayList<Integer>();
	  findFreqPeaks(f, peaks);
	  final int N = peaks.size(); 
	  
	  // compute inharmonicity factor
	  float inh = 0;  // result
	  for(int i=0; i<N; i++){
	    int index = (Integer) peaks.get(i);
	    float fn = f.indexToFreq(index);
	    float ini = Math.abs((fn+1) - (i+1)*f0) / ((i+1) * f0);
	    inh += ini;
	  }
	  
	  return(inh);
	}//end getRightInharmonicity

	/**
	 * This function computes the fundamental frequency of the examined signal.
	 * 
	 *   @param x
	 * 			The time-domain signal to be analyzed
	 * 	 @return f0
	 * 			Fundamental Frequency measured in Hertz [Hz]
	 * 
	 *   @see <a href="http://cnx.org/content/m11714/latest/">Pitch Detection algorithm</a>
	 */
	public float getRightFundFreq(AudioSource x){
	  // compute FFT
	  FFT f = new FFT(x.bufferSize(), x.sampleRate());
	  f.window(FFT.HAMMING);
	  f.forward(x.right);
	  final int ss = f.specSize();  // spectrum size
	  
	  // downsample the spectrum by a factor of 2
	  FFT fd2 = new FFT(x.bufferSize(), x.sampleRate() / 2);
	  fd2.window(FFT.HAMMING);
	  fd2.forward(x.right);
	  
	  // downsample the spectrum by a factor of 3
	  FFT fd3 = new FFT(x.bufferSize(), x.sampleRate() / 3);
	  fd3.window(FFT.HAMMING);
	  fd3.forward(x.right);
	  
	  // multiply the trhee spectra together
	  float[] y = new float[ss];
	  for(int i=0; i<ss-1; i++){
	    y[i] = f.getBand(i) * fd2.getBand(i) * fd3.getBand(i);
	  }
	  
	  // find position of the maximum peak in the resulting spectrum
	  float f0Val = -1;
	  int f0Pos = -1;
	  for(int i=0; i<ss; i++){
	    if(y[i] > f0Val){
	      f0Pos = i;
	      f0Val = y[i];
	    }
	  }
	  
	  // convert the index into frequency
	  float f0 = f.indexToFreq(f0Pos);
	  return(f0);
	}//end getFundFreq
	
	/**
	 * The irregularity of a spectrum is the degree of variation of the successive peaks of the spectrum. It is computed as the sum of the square of the difference in amplitude between adjoining partials.
	 * 
	 *   @param x
	 * 			The time-domain signal to be analyzed
	 * 	 @return irr
	 * 			Spectral Irregularity coefficient
	 */
	public float getRightIrregularity(AudioSource x){
	  // compute FFT
	  FFT f = new FFT(x.bufferSize(), x.sampleRate());
	  f.window(FFT.HAMMING);
	  f.forward(x.right);
	  
	  float irr = 0;  // result
	  float a = 0;
	  float b = 0; 
	  float num = 0;
	  float den = 0;
	  
	  // compute the spectrum peaks
	  ArrayList<Integer> peaks = new ArrayList<Integer>();
	  findFreqPeaks(f, peaks);
	  final int N = peaks.size();
	  
	  // compute the irregularity between adjacent peaks
	  int idx1, idx2;  // peaks indexes
	  for(int i=0; i<N-1; i++){
	    idx1 = (int) peaks.get(i);
	    idx2 = (int) peaks.get(i+1);
	    a = f.getBand(idx1);
	    b = f.getBand(idx2);
	    
	    num += (Math.pow(a - b, 2));
	    den += Math.pow(a, 2);
	  }
	  
	  if(den != 0){
		  irr = num / den;
	  }
	  else{
		  irr = -1;
	  }
	  return(irr);
	}//end getRightIrregularity
	
	
	/**
	 * The kurtosis gives a measure of the flatness 
	 * of the spectrum around its mean value. 
	 * It is computed from the 4th order of the moment.
	 * 
	 *   @param x
	 * 			The time-domain signal to be analyzed
	 * 	 @return kurt
	 * 			Spectral kurtosis coefficient
	 */
	public float getRightKurtosis(AudioSource x){
	  // compute FFT
	  FFT f = new FFT(x.bufferSize(), x.sampleRate());
	  f.window(FFT.HAMMING);
	  f.forward(x.right);
	  int ss = f.specSize();  // spectrum size
	  
	  float kurt = 0;  // result
	  float xi = 0;
	  float fi = 0;
	  
	  // compute centroid and spread
	  float sc = getRightCentroid(x);  // spectrum centroid
	  float sigma = (float) Math.sqrt(getRightSpread(x));  // spectrum standard deviation
	  
	  if(sigma != 0){
		  // compute kurtosis
		  for(int i=0; i<ss-1; i++){  // for each frequency band...
			  xi = f.getBand(i);
			  fi = f.indexToFreq(i);
			  kurt += (xi * (Math.pow((fi - sc), 4)));  // Why 4? By definition...
		  }

		  kurt /= Math.pow(sigma, 4);
	  }
	  else{
		  kurt = -1;
	  }
	  return(kurt);
	}//end getRightKurtosis
	
	
	/**
	 * Mel-Frequency Cepstral Coefficients (MFCC) originate 
	 * from automatic speech recognition.
	 * The MFCCs are to some extent created according to 
	 * the principles of the human auditory system, but also 
	 * to be a compact representation of the amplitude 
	 * spectrum and with considerations of the computational 
	 * complexity.
	 * 
	 * @param x
	 * 			The time-domain signal to be analyzed
	 * @param index
	 * 			The order of the MFCC to be computed (usually 1 ² <code>index</code> ² 13)
	 * @return melCoeff
	 * 			The <code>index</code>-th order Mel Cepstral coefficient
	 * @see <a href="http://www.speech-recognition.de/matlab-files.html">Mel filterbank algorithm</a>
	 */

	public float getRightMFCC(AudioSource x, int index){
	  // compute FFT
	  FFT f = new FFT(x.bufferSize(), x.sampleRate());
	  f.window(FFT.HAMMING);
	  f.forward(x.right);
	  final int N = f.specSize();  // spectrum size
	  
	  float fs = x.sampleRate();  // sampling frequency
	  final int nofChannels = 22;  // number of mel channels
	  
	  // compute resolution etc
	  final float df = fs / N;  //  frequency increment on linear scale
	  final float Nmax = N; //  maximum fft index
	  final float fmax = fs/2; //  maximum frequency
	  final float melmax = freq2mel(fmax); //  maximum mel frequency
	  final float melinc = melmax / (nofChannels + 1); //  frequency increment on mel scale
	  
	  //  compute center frequencies on mel scale
	  float[] melcenters = new float[nofChannels];
	  for(int i=0; i<nofChannels; i++){
	    melcenters[i] = (i + 1) * melinc;
	  }
	  
	  // compute center frequencies in linear scale [Hz]
	  float[] fcenters = new float[nofChannels];
	  for(int i=0; i<nofChannels; i++){
	    fcenters[i] = mel2freq(melcenters[i]);
	  }
	  
	  // compute indexes of center frequencies
	  float[] indexcenter = new float[nofChannels];
	  for(int i=0; i<nofChannels; i++){
	    indexcenter[i] = Math.round(fcenters[i] / df);
	  }
	  
	  float[] indexstart = new float[nofChannels];  // compute start indices of windows
	  float[] indexstop = new float[nofChannels];  // compute stop indices of windows
	  float[] idxbw = new float[nofChannels];  // compute bandwidth (number of indices per window)
	  for(int i=0; i<nofChannels; i++){
	    if(i == 0){
	      indexstart[i] = 0;
	    }
	    else if(i == nofChannels - 1){
	      indexstop[i] = Nmax-1;
	    }
	    else{
	      indexstart[i] = indexcenter[i-1];
	      indexstop[i] = indexcenter[i+1];
	    }
	    idxbw[i] = (indexstop[i] - indexstart[i]) + 1;
	  }//end for
	  
	  // signal spectrum in mel components
	  float[] melSpectrum = new float[N];
	  for(int i=0; i<N; i++){
	    melSpectrum[i] = 0;
	  }
	  for(int i=0; i<N; i++){
	    int currMelChannel = 0;
	    for(int j=0; j<nofChannels-1; j++){
	      if(i >= indexstart[j] && i < indexstart[j+1]){
	        currMelChannel = j;  // current mel channel
	      }
	    }
	    
	    float Si = f.getBand(i);
	    float Nmel = idxbw[currMelChannel];
	    float currMelIndex = indexstop[currMelChannel] - i;
	    melSpectrum[i] += Si * (2.0 / (Nmel - 1.0) * (((Nmel - 1.0) / 2.0) - Math.abs(currMelIndex - ((Nmel - 1.0) / 2.0))));
	  }//end for
	  
	  // take the log of the mel frequency components
	  for(int i=0; i<N; i++){
	    if(melSpectrum[i] > 0){  // avoid -Infinity
	      float temp = (float) Math.log(melSpectrum[i]);
	      melSpectrum[i] = temp;
	    }
	    else{
	      melSpectrum[i] = 0;
	    }
	  }
	  
	  // perform DCT on the log values
	  float melCoeff = 0;
	  for(int i=0; i<N-1; i++){
	    melCoeff += melSpectrum[i] * Math.cos((Math.PI/N) * (i + 0.5) * index);
	  }
	  
	  return(melCoeff);
	}//end getRightMFCC

	/**
	 * Convert a frequency component into mel scale
	 * 
	 * @param fin
	 * 			The input spectrum frequency
	 * @return mel
	 * 			The corresponding mel-scale component
	 */
	private float freq2mel(float fin){
	  float mel = (float) (1127.0148 * Math.log(1.0 + (fin / 700.0)));
	  return(mel);
	}//end freq2mel

	/**
	 * Convert a mel component into a frequency
	 * 
	 * @param mel
	 * 			The input mel-scale component
	 * @return fout
	 * 			The corresponding spectrum frequency
	 */
	private float mel2freq(float mel){
	  float fout = (float) (700.0 * ((Math.exp(mel / 1127.0148)) - 1.0));
	  return(fout);
	}//end mel2freq
	
	
	/**
	 * The spectrum spread is also called instantaneous 
	 * bandwidth; it is defined as the second central 
	 * moment (i.e. the variance) of spectrum.
	 * 
	 *   @param x
	 * 			The time-domain signal to be analyzed
	 * 	 @return ass
	 * 			The Audio Spectrum Spread coefficient
	 */
	public float getRightSpread(AudioSource x){
	  // compute FFT
	  FFT f = new FFT(x.bufferSize(), x.sampleRate());
	  f.window(FFT.HAMMING);
	  f.forward(x.right);
	  final int Nft = f.specSize();  // spectrum size
	  
	  float deltaF = x.sampleRate() / Nft;  // frequency interval between two FFT bins
	  int kLow = (int) Math.floor(62.5 / deltaF);  // preventing disproportionate weight below this frequency
	  
	  float asc = getRightCentroid(x);  // compute spectrum centroid
	  
	  float ass = 0;  // result
	  float num = 0;
	  float den = 0;
	  
	  // compute first component
	  float Pk = 0;
	  if(kLow > 0){
	    for(int i=0; i<kLow; i++){  // for each band
	      Pk += f.getBand(i);
	    }
	    num += (31.25 - asc) * Pk;  // 31.25 Hz = fLow
	    den += Pk;
	  }
	  
	  // compute the remaining components
	  float fk = 0;
	  for(int i=kLow; i<Nft/2; i++){
	    Pk = f.getBand(i);
	    fk = f.indexToFreq(i);
	    num += (fk - asc) * Pk;
	    den += Pk;
	  }
	  
	  ass = num / den;
	  return(ass);
	}//end getRightSpread
	
	
	/**
	 * Short Time Energy (STE) expresses the instant energy of each frame, and is commonly known as volume.
	 * 
	 * A speech audio sample will give large variations in
	 * the STE-values. Music has much shorter pauses, or no 
	 * pauses, and will therefore have a more constant 
	 * STE-level.
	 *   @param x
	 * 			The time-domain signal to be analyzed
	 *   @return ste
	 *   		Short Time Energy coefficient
	 */
	public float getRightSTE(AudioSource x){
	  final int L = x.bufferSize();
	  float ste = 0;  // result
	  float a;
	  
	  // compute STE
	  for(int i=0; i<L; i++){
	    a = x.right.get(i);
	    ste += (a * a);
	  }
	  
	  return(ste);
	  
	}//end getRightSTE
	
	
	
	/**
	 * The zero-crossing rate is the rate of sign-changes
	 * along a signal, i.e., the rate at which the signal 
	 * changes from positive to negative or back.
	 * 
	 *   @param x
	 * 			The time-domain signal to be analyzed
	 * @return z
	 * 			The amount of sign-changes within the signal
	 */
	public int getRightZCR(AudioSource x){
	  final int L = x.bufferSize();
	  int z = 0;  // result
	  float a, b;
	  
	  // compute ZCR
	  for(int i=0; i<L-1; i++){
	    a = x.right.get(i);
	    b = x.right.get(i+1);
	    if(a*b < 0){  // different signs
	      z++;
	    }
	  }
	  
	  return(z);
	}//end getRightZCR
	
	
	/**
	 * The skewness gives a measure of the asymmetry of the 
	 * spectrum around its mean value. It is computed from 
	 * the 3rd order of the moment. 
	 *   @param x
	 * 			The time-domain signal to be analyzed
	 * @return skew
	 * 			Spectral Skewness coefficient
	 */
	public float getRightSkewness(AudioSource x){
	  // compute FFT
	  FFT f = new FFT(x.bufferSize(), x.sampleRate());
	  f.window(FFT.HAMMING);
	  f.forward(x.right);
	  int ss = f.specSize();  // spectrum size
	  
	  float skew = 0;  // result
	  float xi = 0;
	  float fi = 0;
	  
	  // compute spectrum centroid and "standard deviation"
	  float sc = getRightCentroid(x);  // spectrum centroid
	  float sigma = (float) Math.sqrt(getRightSpread(x));  // spectrum standard deviation
	  
	  // compute skewness
	  for(int i=0; i<ss-1; i++){  // for each frequency band...
	    xi = f.getBand(i);
	    fi = f.indexToFreq(i);
	    skew += (xi * (Math.pow((fi - sc), 3)));  // Why 3? By definition...
	  }
	  
	  skew /= Math.pow(sigma, 3);
	  return(skew);
	}//end getRightSkewness
	
	
	/**
	 * The roughness of a signal is related to the 
	 * beating phenomena between the peaks of its 
	 * spectrum.
	 * 
	 * An estimation of the total roughness can be 
	 * obtained by computing the peaks of the 
	 * spectrum, and taking the average of all the 
	 * dissonance between all possible pairs of 
	 * peaks.
	 *
	 *   @param x
	 * 			The time-domain signal to be analyzed
	 * 	 @return r
	 * 			Roughness coefficient
	 */
	public float getRightRoughness(AudioSource x){
	  // compute FFT
	  FFT f = new FFT(x.bufferSize(), x.sampleRate());
	  f.window(FFT.HAMMING);
	  f.forward(x.right);
	  ArrayList<Integer> peaks = new ArrayList<Integer>();  // spectrum peaks
	  
	  float a = 0;
	  float b = 0;
	  int aPos, bPos;
	  float Amin, Amax, fmin, fmax, X, Y, Z;
	  
	  // constant coefficients
	  final float b1 = (float) 3.5;
	  final float b2 = (float) 5.75;
	  final float s1 = (float) 0.0207;
	  final float s2 = (float) 18.96;
	  
	  float r = 0;  // result
	  
	  // find frequency peaks
	  findFreqPeaks(f, peaks);
	  final int K = peaks.size();
	  
	  for(int i=0; i<K-1; i++){
	    for(int j=i+1; j<K; j++){
	      aPos = (Integer) peaks.get(i);
	      bPos = (Integer) peaks.get(j);
	      a = f.getBand(aPos);
	      b = f.getBand(bPos);
	      Amin = Math.min(a, b);  // minimum frequency component
	      Amax = Math.max(a, b);  // maximum frequency component
	      fmin = f.indexToFreq(aPos);
	      fmax = f.indexToFreq(bPos);
	      
	      final float s = (float) (0.24 / (s1 * fmin + s2));
	      
	      X = Amin * Amax;
	      Y = (2 * Amin) / (Amin + Amax);
	      Z = (float) (Math.exp(-b1 * s * (fmax - fmin)) - Math.exp(-b2 * s * (fmax - fmin)));
	      
	      r += (Math.pow(X, 0.1) * 0.5 * (Math.pow(Y, 3.11)) * Z);
	    }//end for
	  }//end for
	  
	  r /= K;
	  return(r);
	}//end getRightRoughness

	
	/**
	 * The spectral brightness is defined as the 
	 * amount of energy in the spectrum starting 
	 * from a specified cutoff frequency 
	 * (in this case 1500 Hz), divided by the 
	 * total energy of the signal spectrum
	 * 
	 *   @param x
	 * 			The time-domain signal to be analyzed
	 * 	 @return b
	 * 			Spectral Brightness coefficient
	 */
	public float getLeftBrightness(AudioSource x){
	  // compute FFT
	  FFT f = new FFT(x.bufferSize(), x.sampleRate());
	  f.window(FFT.HAMMING);
	  f.forward(x.left);
	  int ss = f.specSize();  // spectrum size
	  
	  final int fc = 1500;  // Why 1500 Hz? By definition...
	  final int kc = f.freqToIndex(fc);
	  
	  float num = 0;  // numerator
	  float den = 0;  // denominator
	  float b = 0;    // result
	  
	  for(int i=0; i<ss; i++){
	    den += f.getBand(i);
	    if(i >= kc){
	      num += f.getBand(i);
	    }
	  }
	  
	  b = num / den;
	  return(b);
	}//end getLeftBrightness


	/**
	 * The spectral centroid is the baricenter 
	 * of the spectrum; it is calculated as the 
	 * weighted mean of the frequencies present 
	 * in the signal, with their magnitude 
	 * as the weights.
	 * 
	 * The result is expressed in Hertz [Hz]
	 * 
	 *   @param x
	 * 			The time-domain signal to be analyzed
	 * 	 @return asc
	 * 			The Audio Spectral Centroid measured in Hertz [Hz]
	 */
	public float getLeftCentroid(AudioSource x){
	  // compute FFT
	  FFT f = new FFT(x.bufferSize(), x.sampleRate());
	  f.window(FFT.HAMMING);
	  f.forward(x.left);
	  int Nft = f.specSize();  // spectrum size
	  
	  float deltaF = x.sampleRate() / Nft;  // frequency interval between two FFT bins
	  int kLow = (int) Math.floor(62.5 / deltaF);  // preventing disproportionate weight below this frequency
	  
	  float asc = 0;  // result
	  float num = 0;
	  float den = 0;
	  
	  // compute the first component
	  float Pk = 0;
	  if(kLow > 0){
	    for(int i=0; i<kLow; i++){  // for each band
	      Pk += f.getBand(i);
	    }
	    num += 31.25 * Pk;  // 31.25 Hz = fLow
	    den += Pk;
	  }
	  
	  // compute the remaining components
	  float fk = 0;
	  for(int i=kLow; i<Nft/2; i++){
	    Pk = f.getBand(i);
	    fk = f.indexToFreq(i);
	    num += fk * Pk;
	    den += Pk;
	  }
	  
	  asc = num / den;
	  return(asc);
	  
	}//end getLeftCentroid


	/**
	 * The Spectral Entropy is computed as the 
	 * ShannonÕs entropy of the signal spectrum.
	 * 
	 * The entropy is divided by the logarithm of 
	 * the length N of the spectrum in order to have 
	 * a result which is independent from the windowing length.
	 * 
	 *   @param x
	 * 			The time-domain signal to be analyzed
	 *	 @return ent
	 *			Spectral Entropy coefficient
	 */
	public float getLeftEntropy(AudioSource x){
	  // compute FFT
	  FFT f = new FFT(x.bufferSize(), x.sampleRate());
	  f.window(FFT.HAMMING);
	  f.forward(x.left);
	  int ss = f.specSize();  // spectrum size
	  
	  float ent = 0;  // result
	  float num = 0;
	  float den = 0;
	  
	  float Xi = 0;
	  for(int i=0; i<ss; i++){
	    Xi = f.getBand(i);
	    num += (Xi * Math.log(Xi));
	  }
	  
	  den = (float) Math.log(ss);  // logarithm of the spectrum length
	  
	  ent = - (num / den);
	  return(ent);
	}//end getLeftEntropy


	/**
	 * Spectral flatness provides a way to quantify 
	 * how tone-like a sound is, as opposed to 
	 * being noise-like. The meaning of tonal in this 
	 * context is in the sense of the amount of peaks 
	 * or resonant structure in a power spectrum, as 
	 * opposed to flat spectrum of a white noise.
	 * 
	 *   @param x
	 * 			The time-domain signal to be analyzed
	 * 	 @return asf
	 * 			Audio Spectral Flatness coefficient
	 */
	public float getLeftFlatness(AudioSource x){
	  // compute FFT
	  FFT f = new FFT(x.bufferSize(), x.sampleRate());
	  f.window(FFT.HAMMING);
	  f.forward(x.left);
	  
	  float num = 1;
	  float den = 0;
	  float Si = 0;
	  float asf = 0;  // result

	  final int n = -8;
	  final int B = 24;  // number of bands
	  final float loF = (float) (Math.pow(2, n/4.0) * 1000);  // lowest frequency [Hz]
	  final float hiF = (float) (Math.pow(2, B/4.0) * loF);  // highest frequency [Hz]
	  final int loK = f.freqToIndex(loF);
	  final int hiK = f.freqToIndex(hiF);
	  final float reduceFactor = hiK - loK + 1;
	  
	  for(int k=loK; k<hiK; k++){
	    Si = f.getBand(k);
	    num *= Math.pow(Si, 1.0/reduceFactor);  // n-th root
	    den += Si;
	  }
	  
	  den /= reduceFactor;
	  asf = num / den;
	  return(asf);
	}//end getLeftFlatness


	/**
	 * Spectral flux is a measure of how quickly 
	 * the power spectrum of a signal is changing, 
	 * calculated by comparing the power spectrum 
	 * of half a frame against the power spectrum 
	 * of the other half.
	 * 
	 *   @param x
	 * 			The time-domain signal to be analyzed
	 *	 @return flux
	 *			Spectral Flux coefficient
	 */
	public float getLeftFlux(AudioSource x){
	  final int bufSize = x.bufferSize();
	  float[] ab = new float[bufSize];
	  float[] ab1, ab2;
	  ab1 = new float[bufSize / 2];
	  ab2 = new float[bufSize / 2];
	  
	  // split the audio buffer into 2 halves
	  ab = x.left.toArray();
	  for(int i=0; i<(bufSize / 2); i++){
	    ab1[i] = ab[i];
	    ab2[i] = ab[(bufSize / 2) + i];
	  }
	  
	  // compute the 2 spectra
	  FFT f1 = new FFT(ab1.length, x.sampleRate());
	  FFT f2 = new FFT(ab2.length, x.sampleRate());
	  f1.window(FFT.HAMMING);
	  f1.forward(ab1);
	  f2.window(FFT.HAMMING);
	  f2. forward(ab2);
	  int ss = f1.specSize();  // spectrum size
	  
	  // compute the spectral flux
	  float S1 = 0;
	  float S2 = 0;
	  for(int i=0; i<ss; i++){
	    S1 += Math.pow(f1.getBand(i), 2);
	    S2 += Math.pow(f2.getBand(i), 2);
	  }
	  
	  float flux = (float) Math.sqrt(Math.pow(S2 - S1, 2));
	  return(flux);
	}//end getLeftFlux


	/**
	 * The harmonic ratio (HR) is a measure of the proportion
	 * of harmonic components in the power spectrum. It is in
	 * practice a coefficient associated to each signal frame,
	 * computed as the maximum value of the autocorrelation 
	 * function.
	 * 
	 *   @param x
	 * 			The time-domain signal to be analyzed
	 * 	 @return hr
	 * 			Harmonic Ratio coefficient
	 */
	public float getLeftHR(AudioSource x){
	  final int L = x.bufferSize();
	  final int Mmin = 5;  // prevent giving higher importance to low lags
	  final int Mmax = L -1;  // maximum lag
	  
	  float hr = -1;  // result
	  float gamma;  // current autocorrelation value
	  // compute autocorrelation
	  for(int m=Mmin; m<=Mmax; m++){
	    gamma = leftAutoCorr(m, x);
	    if(gamma > hr){
	      hr = gamma;
	    }
	  }
	  
	  return(hr);
	}//end getLeftZCR

	/**
	 * Autocorrelation function
	 * 
	 *   @param m
	 *   		Imput lag
	 *   @param x
	 * 			The time-domain signal to be analyzed
	 * @return gamma
	 * 			<code>m</code>-lag autocorrelation output
	 * 
	 * @see <a href="http://en.wikipedia.org/wiki/Autocorrelation">Definition of autocorrelation</a>
	 */
	private float leftAutoCorr(int m, AudioSource x){
	  final int L = x.bufferSize();
	  float num = 0;
	  float den1 = 0;
	  float den2 = 0;
	  float gamma = 0;  // result
	  
	  // compute autocorrelation
	  for(int t=0; t<L-1; t++){
	    if(t+m < L){
	      num += x.left.get(t) * x.left.get(t+m);
	      den1 += Math.pow(x.left.get(t), 2);
	      den2 += Math.pow(x.left.get(t+m), 2);
	    }
	  }
	  
	  gamma = (float) (num / Math.sqrt(den1 * den2));
	  return(gamma);
	}//end autoCorr


	/**
	 * Inharmonicity estimates the amount of partials that 
	 * are not multiples of the fundamental frequency. It is 
	 * computed by considering the ideal locations of 
	 * harmonics vs the actual harmonics in a spectrum.
	 * 
	 * In this implementation, a lower value implies a more
	 * harmonic signal.
	 * 
	 *   @param x
	 * 			The time-domain signal to be analyzed
	 * 	 @return inh
	 * 			Inharmonicity coefficient
	 */
	public float getLeftInharmonicity(AudioSource x){
	  // compute FFT
	  FFT f = new FFT(x.bufferSize(), x.sampleRate());
	  f.window(FFT.HAMMING);
	  f.forward(x.left);
	  
	  // compute the fundamental frequency f0
	  float f0 = getLeftFundFreq(x);
	  
	  // compute the spectrum peaks
	  ArrayList<Integer> peaks = new ArrayList<Integer>();
	  findFreqPeaks(f, peaks);
	  final int N = peaks.size(); 
	  
	  // compute inharmonicity factor
	  float inh = 0;  // result
	  for(int i=0; i<N; i++){
	    int index = (Integer) peaks.get(i);
	    float fn = f.indexToFreq(index);
	    float ini = Math.abs((fn+1) - (i+1)*f0) / ((i+1) * f0);
	    inh += ini;
	  }
	  
	  return(inh);
	}//end getLeftInharmonicity

	/**
	 * This function computes the fundamental frequency of the examined signal.
	 *   @param x
	 * 			The time-domain signal to be analyzed
	 * 	 @return f0
	 * 			Fundamental Frequency measured in Hertz [Hz]
	 * 
	 *   @see <a href="http://cnx.org/content/m11714/latest/">Pitch Detection algorithm</a>
	 */
	public float getLeftFundFreq(AudioSource x){
	  // compute FFT
	  FFT f = new FFT(x.bufferSize(), x.sampleRate());
	  f.window(FFT.HAMMING);
	  f.forward(x.left);
	  final int ss = f.specSize();  // spectrum size
	  
	  // downsample the spectrum by a factor of 2
	  FFT fd2 = new FFT(x.bufferSize(), x.sampleRate() / 2);
	  fd2.window(FFT.HAMMING);
	  fd2.forward(x.left);
	  
	  // downsample the spectrum by a factor of 3
	  FFT fd3 = new FFT(x.bufferSize(), x.sampleRate() / 3);
	  fd3.window(FFT.HAMMING);
	  fd3.forward(x.left);
	  
	  // multiply the trhee spectra together
	  float[] y = new float[ss];
	  for(int i=0; i<ss-1; i++){
	    y[i] = f.getBand(i) * fd2.getBand(i) * fd3.getBand(i);
	  }
	  
	  // find position of the maximum peak in the resulting spectrum
	  float f0Val = -1;
	  int f0Pos = -1;
	  for(int i=0; i<ss; i++){
	    if(y[i] > f0Val){
	      f0Pos = i;
	      f0Val = y[i];
	    }
	  }
	  
	  // convert the index into frequency
	  float f0 = f.indexToFreq(f0Pos);
	  return(f0);
	}//end getFundFreq

	/**
	 * The irregularity of a spectrum is the degree of variation of the successive peaks of the spectrum. It is computed as the sum of the square of the difference in amplitude between adjoining partials.
	 * 
	 *   @param x
	 * 			The time-domain signal to be analyzed
	 *   @return irr
	 *   		Spectrum Irregularity coefficient
	 */
	public float getLeftIrregularity(AudioSource x){
	  // compute FFT
	  FFT f = new FFT(x.bufferSize(), x.sampleRate());
	  f.window(FFT.HAMMING);
	  f.forward(x.left);
	  
	  float irr = 0;  // result
	  float a = 0;
	  float b = 0; 
	  float num = 0;
	  float den = 0;
	  
	  // compute the spectrum peaks
	  ArrayList<Integer> peaks = new ArrayList<Integer>();
	  findFreqPeaks(f, peaks);
	  final int N = peaks.size();
	  
	  // compute the irregularity between adjacent peaks
	  int idx1, idx2;  // peaks indexes
	  for(int i=0; i<N-1; i++){
	    idx1 = (int) peaks.get(i);
	    idx2 = (int) peaks.get(i+1);
	    a = f.getBand(idx1);
	    b = f.getBand(idx2);
	    
	    num += (Math.pow(a - b, 2));
	    den += Math.pow(a, 2);
	  }
	  
	  if(den != 0){
		  irr = num / den;
	  }
	  else{
		  irr = -1;
	  }
	  return(irr);
	}//end getLeftIrregularity


	/**
	 * The kurtosis gives a measure of the flatness 
	 *  of the spectrum around its mean value. 
	 *  It is computed from the 4th order of the moment.
	 * 
	 *   @param x
	 * 			The time-domain signal to be analyzed
	 *   @return kurt
	 *   		Spectral Kurtosis coefficient
	 */
	public float getLeftKurtosis(AudioSource x){
	  // compute FFT
	  FFT f = new FFT(x.bufferSize(), x.sampleRate());
	  f.window(FFT.HAMMING);
	  f.forward(x.left);
	  int ss = f.specSize();  // spectrum size
	  
	  float kurt = 0;  // result
	  float xi = 0;
	  float fi = 0;
	  
	  // compute centroid and spread
	  float sc = getLeftCentroid(x);  // spectrum centroid
	  float sigma = (float) Math.sqrt(getLeftSpread(x));  // spectrum standard deviation
	  
	  if(sigma != 0){
		  // compute kurtosis
		  for(int i=0; i<ss-1; i++){  // for each frequency band...
			  xi = f.getBand(i);
			  fi = f.indexToFreq(i);
			  kurt += (xi * (Math.pow((fi - sc), 4)));  // Why 4? By definition...
		  }

		  kurt /= Math.pow(sigma, 4);
	  }
	  else{
		  kurt = -1;
	  }
	  return(kurt);
	}//end getLeftKurtosis
	
	
	/**
	 * Mel-Frequency Cepstral Coefficients (MFCC) originate 
	 * from automatic speech recognition.
	 * The MFCCs are to some extent created according to 
	 * the principles of the human auditory system, but also 
	 * to be a compact representation of the amplitude 
	 * spectrum and with considerations of the computational 
	 * complexity.
	 * 
	 * @param x
	 * 			The time-domain signal to be analyzed
	 * @param index
	 * 			The order of the MFCC to be computed (usually 1 ² <code>index</code> ² 13)
	 * @return melCoeff
	 * 			The <code>index</code>-th order Mel Cepstral coefficient
	 * @see <a href="http://www.speech-recognition.de/matlab-files.html">Mel filterbank algorithm</a>
	 */
	public float getLeftMFCC(AudioSource x, int index){
	  // compute FFT
	  FFT f = new FFT(x.bufferSize(), x.sampleRate());
	  f.window(FFT.HAMMING);
	  f.forward(x.left);
	  final int N = f.specSize();  // spectrum size
	  
	  float fs = x.sampleRate();  // sampling frequency
	  final int nofChannels = 22;  // number of mel channels
	  
	  // compute resolution etc
	  final float df = fs / N;  //  frequency increment on linear scale
	  final float Nmax = N; //  maximum fft index
	  final float fmax = fs/2; //  maximum frequency
	  final float melmax = freq2mel(fmax); //  maximum mel frequency
	  final float melinc = melmax / (nofChannels + 1); //  frequency increment on mel scale
	  
	  //  compute center frequencies on mel scale
	  float[] melcenters = new float[nofChannels];
	  for(int i=0; i<nofChannels; i++){
	    melcenters[i] = (i + 1) * melinc;
	  }
	  
	  // compute center frequencies in linear scale [Hz]
	  float[] fcenters = new float[nofChannels];
	  for(int i=0; i<nofChannels; i++){
	    fcenters[i] = mel2freq(melcenters[i]);
	  }
	  
	  // compute indexes of center frequencies
	  float[] indexcenter = new float[nofChannels];
	  for(int i=0; i<nofChannels; i++){
	    indexcenter[i] = Math.round(fcenters[i] / df);
	  }
	  
	  float[] indexstart = new float[nofChannels];  // compute start indices of windows
	  float[] indexstop = new float[nofChannels];  // compute stop indices of windows
	  float[] idxbw = new float[nofChannels];  // compute bandwidth (number of indices per window)
	  for(int i=0; i<nofChannels; i++){
	    if(i == 0){
	      indexstart[i] = 0;
	    }
	    else if(i == nofChannels - 1){
	      indexstop[i] = Nmax-1;
	    }
	    else{
	      indexstart[i] = indexcenter[i-1];
	      indexstop[i] = indexcenter[i+1];
	    }
	    idxbw[i] = (indexstop[i] - indexstart[i]) + 1;
	  }//end for
	  
	  // signal spectrum in mel components
	  float[] melSpectrum = new float[N];
	  for(int i=0; i<N; i++){
	    melSpectrum[i] = 0;
	  }
	  for(int i=0; i<N; i++){
	    int currMelChannel = 0;
	    for(int j=0; j<nofChannels-1; j++){
	      if(i >= indexstart[j] && i < indexstart[j+1]){
	        currMelChannel = j;  // current mel channel
	      }
	    }
	    
	    float Si = f.getBand(i);
	    float Nmel = idxbw[currMelChannel];
	    float currMelIndex = indexstop[currMelChannel] - i;
	    melSpectrum[i] += Si * (2.0 / (Nmel - 1.0) * (((Nmel - 1.0) / 2.0) - Math.abs(currMelIndex - ((Nmel - 1.0) / 2.0))));
	  }//end for
	  
	  // take the log of the mel frequency components
	  for(int i=0; i<N; i++){
	    if(melSpectrum[i] > 0){  // avoid -Infinity
	      float temp = (float) Math.log(melSpectrum[i]);
	      melSpectrum[i] = temp;
	    }
	    else{
	      melSpectrum[i] = 0;
	    }
	  }
	  
	  // perform DCT on the log values
	  float melCoeff = 0;
	  for(int i=0; i<N-1; i++){
	    melCoeff += melSpectrum[i] * Math.cos((Math.PI/N) * (i + 0.5) * index);
	  }
	  
	  return(melCoeff);
	}//end getLeftMFCC


	/**
	 * The spectrum spread is also called instantaneous 
	 *  bandwidth; it is defined as the second central 
	 *  moment (i.e. the variance) of spectrum.
	 * 
	 *   @param x
	 * 			The time-domain signal to be analyzed
	 * 	 @return ass
	 * 			Audio Spectrum Spread coefficient
	 */

	public float getLeftSpread(AudioSource x){
	  // compute FFT
	  FFT f = new FFT(x.bufferSize(), x.sampleRate());
	  f.window(FFT.HAMMING);
	  f.forward(x.left);
	  final int Nft = f.specSize();  // spectrum size
	  
	  float deltaF = x.sampleRate() / Nft;  // frequency interval between two FFT bins
	  int kLow = (int) Math.floor(62.5 / deltaF);  // preventing disproportionate weight below this frequency
	  
	  float asc = getLeftCentroid(x);  // compute spectrum centroid
	  
	  float ass = 0;  // result
	  float num = 0;
	  float den = 0;
	  
	  // compute first component
	  float Pk = 0;
	  if(kLow > 0){
	    for(int i=0; i<kLow; i++){  // for each band
	      Pk += f.getBand(i);
	    }
	    num += (31.25 - asc) * Pk;  // 31.25 Hz = fLow
	    den += Pk;
	  }
	  
	  // compute the remaining components
	  float fk = 0;
	  for(int i=kLow; i<Nft/2; i++){
	    Pk = f.getBand(i);
	    fk = f.indexToFreq(i);
	    num += (fk - asc) * Pk;
	    den += Pk;
	  }
	  
	  ass = num / den;
	  return(ass);
	}//end getLeftSpread


	/**
	 * Short Time Energy (STE) expresses the instant energy of each frame, and is commonly known as volume.
	 * 
	 * A speech audio sample will give large variations in
	 * the STE-values. Music has much shorter pauses, or no 
	 * pauses, and will therefore have a more constant 
	 * STE-level.
	 *   @param x
	 * 			The time-domain signal to be analyzed
	 *   @return ste
	 *   		Short Time Energy coefficient
	 */
	public float getLeftSTE(AudioSource x){
	  final int L = x.bufferSize();
	  float ste = 0;  // result
	  float a;
	  
	  // compute STE
	  for(int i=0; i<L; i++){
	    a = x.left.get(i);
	    ste += (a * a);
	  }
	  
	  return(ste);
	  
	}//end getLeftSTE



	/**
	 * The zero-crossing rate is the rate of sign-changes
	 * along a signal, i.e., the rate at which the signal 
	 * changes from positive to negative or back.
	 * 
	 *   @param x
	 * 			The time-domain signal to be analyzed
	 *	 @return z
	 *			The amount of sign changes within the input signal
	 */

	public int getLeftZCR(AudioSource x){
	  final int L = x.bufferSize();
	  int z = 0;  // result
	  float a, b;
	  
	  // compute ZCR
	  for(int i=0; i<L-1; i++){
	    a = x.left.get(i);
	    b = x.left.get(i+1);
	    if(a*b < 0){  // different signs
	      z++;
	    }
	  }
	  
	  return(z);
	}//end getLeftZCR


	/**
	 * The skewness gives a measure of the asymmetry of the 
	 * spectrum around its mean value. It is computed from 
	 * the 3rd order of the moment.
	 *  
	 *   @param x
	 * 			The time-domain signal to be analyzed
	 *	 @return skew
	 *			Spectrum Skewness coefficient
	 */

	public float getLeftSkewness(AudioSource x){
	  // compute FFT
	  FFT f = new FFT(x.bufferSize(), x.sampleRate());
	  f.window(FFT.HAMMING);
	  f.forward(x.left);
	  int ss = f.specSize();  // spectrum size
	  
	  float skew = 0;  // result
	  float xi = 0;
	  float fi = 0;
	  
	  // compute spectrum centroid and "standard deviation"
	  float sc = getLeftCentroid(x);  // spectrum centroid
	  float sigma = (float) Math.sqrt(getLeftSpread(x));  // spectrum standard deviation
	  
	  // compute skewness
	  for(int i=0; i<ss-1; i++){  // for each frequency band...
	    xi = f.getBand(i);
	    fi = f.indexToFreq(i);
	    skew += (xi * (Math.pow((fi - sc), 3)));  // Why 3? By definition...
	  }
	  
	  skew /= Math.pow(sigma, 3);
	  return(skew);
	}//end getLeftSkewness


	/**
	 * The roughness of a signal is related to the 
	 * beating phenomena between the peaks of its 
	 * spectrum.
	 * 
	 * An estimation of the total roughness can be 
	 * obtained by computing the peaks of the 
	 * spectrum, and taking the average of all the 
	 * dissonance between all possible pairs of 
	 * peaks.
	 * 
	 *   @param x
	 * 			The time-domain signal to be analyzed
	 *	 @return r
	 *			Roughness coefficient
	 */

	public float getLeftRoughness(AudioSource x){
	  // compute FFT
	  FFT f = new FFT(x.bufferSize(), x.sampleRate());
	  f.window(FFT.HAMMING);
	  f.forward(x.left);
	  ArrayList<Integer> peaks = new ArrayList<Integer>();  // spectrum peaks
	  
	  float a = 0;
	  float b = 0;
	  int aPos, bPos;
	  float Amin, Amax, fmin, fmax, X, Y, Z;
	  
	  // constant coefficients
	  final float b1 = (float) 3.5;
	  final float b2 = (float) 5.75;
	  final float s1 = (float) 0.0207;
	  final float s2 = (float) 18.96;
	  
	  float r = 0;  // result
	  
	  // find frequency peaks
	  findFreqPeaks(f, peaks);
	  final int K = peaks.size();
	  
	  for(int i=0; i<K-1; i++){
	    for(int j=i+1; j<K; j++){
	      aPos = (Integer) peaks.get(i);
	      bPos = (Integer) peaks.get(j);
	      a = f.getBand(aPos);
	      b = f.getBand(bPos);
	      Amin = Math.min(a, b);  // minimum frequency component
	      Amax = Math.max(a, b);  // maximum frequency component
	      fmin = f.indexToFreq(aPos);
	      fmax = f.indexToFreq(bPos);
	      
	      final float s = (float) (0.24 / (s1 * fmin + s2));
	      
	      X = Amin * Amax;
	      Y = (2 * Amin) / (Amin + Amax);
	      Z = (float) (Math.exp(-b1 * s * (fmax - fmin)) - Math.exp(-b2 * s * (fmax - fmin)));
	      
	      r += (Math.pow(X, 0.1) * 0.5 * (Math.pow(Y, 3.11)) * Z);
	    }//end for
	  }//end for
	  
	  r /= K;
	  return(r);
	}//end getLeftRoughness
	
	
	/**
	 * The spectral brightness is defined as the 
	 * amount of energy in the spectrum starting 
	 * from a specified cutoff frequency 
	 * (in this case 1500 Hz), divided by the 
	 * total energy of the signal spectrum
	 * 
	 *   @param x
	 * 			The time-domain signal to be analyzed
	 * 	 @return b
	 * 			Spectral Brightness coefficient
	 */

	public float getBrightness(AudioSource x){
	  // compute FFT
	  FFT f = new FFT(x.bufferSize(), x.sampleRate());
	  f.window(FFT.HAMMING);
	  f.forward(x.mix);
	  int ss = f.specSize();  // spectrum size
	  
	  final int fc = 1500;  // Why 1500 Hz? By definition...
	  final int kc = f.freqToIndex(fc);
	  
	  float num = 0;  // numerator
	  float den = 0;  // denominator
	  float b = 0;    // result
	  
	  for(int i=0; i<ss; i++){
	    den += f.getBand(i);
	    if(i >= kc){
	      num += f.getBand(i);
	    }
	  }
	  
	  b = num / den;
	  return(b);
	}//end getBrightness
	
	
	/**
	 * The spectral centroid is the baricenter 
	 * of the spectrum; it is calculated as the 
	 * weighted mean of the frequencies present 
	 * in the signal, with their magnitude 
	 * as the weights.
	 * 
	 * The result is expressed in Hertz [Hz]
	 * 
	 *   @param x
	 * 			The time-domain signal to be analyzed
	 *	 @return asc
	 *			Audio Spectral Centroid coefficient
	 */
	public float getCentroid(AudioSource x){
	  // compute FFT
	  FFT f = new FFT(x.bufferSize(), x.sampleRate());
	  f.window(FFT.HAMMING);
	  f.forward(x.mix);
	  int Nft = f.specSize();  // spectrum size
	  
	  float deltaF = x.sampleRate() / Nft;  // frequency interval between two FFT bins
	  int kLow = (int) Math.floor(62.5 / deltaF);  // preventing disproportionate weight below this frequency
	  
	  float asc = 0;  // result
	  float num = 0;
	  float den = 0;
	  
	  // compute the first component
	  float Pk = 0;
	  if(kLow > 0){
	    for(int i=0; i<kLow; i++){  // for each band
	      Pk += f.getBand(i);
	    }
	    num += 31.25 * Pk;  // 31.25 Hz = fLow
	    den += Pk;
	  }
	  
	  // compute the remaining components
	  float fk = 0;
	  for(int i=kLow; i<Nft/2; i++){
	    Pk = f.getBand(i);
	    fk = f.indexToFreq(i);
	    num += fk * Pk;
	    den += Pk;
	  }
	  
	  asc = num / den;
	  return(asc);
	  
	}//end getCentroid
	
	
	/**
	 * The Spectral Entropy is computed as the 
	 * ShannonÕs entropy of the signal spectrum.
	 * 
	 * The entropy is divided by the logarithm of 
	 * the length N of the spectrum in order to have 
	 * a result which is independent from the windowing length.
	 * 
	 *   @param x
	 * 			The time-domain signal to be analyzed
	 *	 @return ent
	 *			Spectral Entropy coefficient
	 */
	public float getEntropy(AudioSource x){
	  // compute FFT
	  FFT f = new FFT(x.bufferSize(), x.sampleRate());
	  f.window(FFT.HAMMING);
	  f.forward(x.mix);
	  int ss = f.specSize();  // spectrum size
	  
	  float ent = 0;  // result
	  float num = 0;
	  float den = 0;
	  
	  float Xi = 0;
	  for(int i=0; i<ss; i++){
	    Xi = f.getBand(i);
	    num += (Xi * Math.log(Xi));
	  }
	  
	  den = (float) Math.log(ss);  // logarithm of the spectrum length
	  
	  ent = - (num / den);
	  return(ent);
	}//end getEntropy
	
	
	/**
	 * Spectral flatness provides a way to quantify 
	 * how tone-like a sound is, as opposed to 
	 * being noise-like. The meaning of tonal in this 
	 * context is in the sense of the amount of peaks 
	 * or resonant structure in a power spectrum, as 
	 * opposed to flat spectrum of a white noise.
	 * 
	 *  @param x
	 * 			The time-domain signal to be analyzed
	 *	 @return asf
	 *			Audio Spectral Flatness coefficient
	 */
	public float getFlatness(AudioSource x){
	  // compute FFT
	  FFT f = new FFT(x.bufferSize(), x.sampleRate());
	  f.window(FFT.HAMMING);
	  f.forward(x.mix);
	  
	  float num = 1;
	  float den = 0;
	  float Si = 0;
	  float asf = 0;  // result

	  final int n = -8;
	  final int B = 24;  // number of bands
	  final float loF = (float) (Math.pow(2, n/4.0) * 1000);  // lowest frequency [Hz]
	  final float hiF = (float) (Math.pow(2, B/4.0) * loF);  // highest frequency [Hz]
	  final int loK = f.freqToIndex(loF);
	  final int hiK = f.freqToIndex(hiF);
	  final float reduceFactor = hiK - loK + 1;
	  
	  for(int k=loK; k<hiK; k++){
	    Si = f.getBand(k);
	    num *= Math.pow(Si, 1.0/reduceFactor);  // n-th root
	    den += Si;
	  }
	  
	  den /= reduceFactor;
	  asf = num / den;
	  return(asf);
	}//end getFlatness
	
	
	/**
	 * Spectral flux is a measure of how quickly 
	 * the power spectrum of a signal is changing, 
	 * calculated by comparing the power spectrum 
	 * of half a frame against the power spectrum 
	 * of the other half.
	 * 
	 *  @param x
	 * 			The time-domain signal to be analyzed
	 * @return flux
	 * 			Spectral Flux coefficient
	 */
	public float getFlux(AudioSource x){
	  final int bufSize = x.bufferSize();
	  float[] ab = new float[bufSize];
	  float[] ab1, ab2;
	  ab1 = new float[bufSize / 2];
	  ab2 = new float[bufSize / 2];
	  
	  // split the audio buffer into 2 halves
	  ab = x.mix.toArray();
	  for(int i=0; i<(bufSize / 2); i++){
	    ab1[i] = ab[i];
	    ab2[i] = ab[(bufSize / 2) + i];
	  }
	  
	  // compute the 2 spectra
	  FFT f1 = new FFT(ab1.length, x.sampleRate());
	  FFT f2 = new FFT(ab2.length, x.sampleRate());
	  f1.window(FFT.HAMMING);
	  f1.forward(ab1);
	  f2.window(FFT.HAMMING);
	  f2. forward(ab2);
	  int ss = f1.specSize();  // spectrum size
	  
	  // compute the spectral flux
	  float S1 = 0;
	  float S2 = 0;
	  for(int i=0; i<ss; i++){
	    S1 += Math.pow(f1.getBand(i), 2);
	    S2 += Math.pow(f2.getBand(i), 2);
	  }
	  
	  float flux = (float) Math.sqrt(Math.pow(S2 - S1, 2));
	  return(flux);
	}//end getFlux
	
	
	/**
	 * The harmonic ratio (HR) is a measure of the proportion
	 * of harmonic components in the power spectrum. It is in
	 * practice a coefficient associated to each signal frame,
	 * computed as the maximum value of the autocorrelation 
	 * function.
	 * 
	 *  @param x
	 * 			The time-domain signal to be analyzed
	 *	 @return hr
	 *			Harmonic Ratio coefficient
	 */
	public float getHR(AudioSource x){
	  final int L = x.bufferSize();
	  final int Mmin = 5;  // prevent giving higher importance to low lags
	  final int Mmax = L -1;  // maximum lag
	  
	  float hr = -1;  // result
	  float gamma;  // current autocorrelation value
	  // compute autocorrelation
	  for(int m=Mmin; m<=Mmax; m++){
	    gamma = autoCorr(m, x);
	    if(gamma > hr){
	      hr = gamma;
	    }
	  }
	  
	  return(hr);
	}//end getZCR

	/**
	 * Autocorrelation function
	 * 
	 *   @param m
	 *   		Imput lag
	 *   @param x
	 * 			The time-domain signal to be analyzed
	 * @return gamma
	 * 			<code>m</code>-lag autocorrelation output
	 * 
	 * @see <a href="http://en.wikipedia.org/wiki/Autocorrelation">Definition of autocorrelation</a>
	 */
	private float autoCorr(int m, AudioSource x){
	  final int L = x.bufferSize();
	  float num = 0;
	  float den1 = 0;
	  float den2 = 0;
	  float gamma = 0;  // result
	  
	  // compute autocorrelation
	  for(int t=0; t<L-1; t++){
	    if(t+m < L){
	      num += x.mix.get(t) * x.mix.get(t+m);
	      den1 += Math.pow(x.mix.get(t), 2);
	      den2 += Math.pow(x.mix.get(t+m), 2);
	    }
	  }
	  
	  gamma = (float) (num / Math.sqrt(den1 * den2));
	  return(gamma);
	}//end autoCorr
	
	
	/**
	 * Inharmonicity estimates the amount of partials that 
	 * are not multiples of the fundamental frequency. It is 
	 * computed by considering the ideal locations of 
	 * harmonics vs the actual harmonics in a spectrum.
	 * 
	 * In this implementation, a lower value implies a more
	 * harmonic signal.
	 * 
	 *  @param x
	 * 			The time-domain signal to be analyzed
	 *  @return inh
	 *  		Inharmonicity coefficient
	 */
	public float getInharmonicity(AudioSource x){
	  // compute FFT
	  FFT f = new FFT(x.bufferSize(), x.sampleRate());
	  f.window(FFT.HAMMING);
	  f.forward(x.mix);
	  
	  // compute the fundamental frequency f0
	  float f0 = getFundFreq(x);
	  
	  // compute the spectrum peaks
	  ArrayList<Integer> peaks = new ArrayList<Integer>();
	  findFreqPeaks(f, peaks);
	  final int N = peaks.size(); 
	  
	  // compute inharmonicity factor
	  float inh = 0;  // result
	  for(int i=0; i<N; i++){
	    int index = (Integer) peaks.get(i);
	    float fn = f.indexToFreq(index);
	    float ini = Math.abs((fn+1) - (i+1)*f0) / ((i+1) * f0);
	    inh += ini;
	  }
	  
	  return(inh);
	}//end getInharmonicity

	/**
	 * This function computes the fundamental frequency of the examined signal.
	 * 
	 *   @param x
	 * 			The time-domain signal to be analyzed
	 * 	 @return f0
	 * 			Fundamental Frequency measured in Hertz [Hz]
	 * 
	 *   @see <a href="http://cnx.org/content/m11714/latest/">Pitch Detection algorithm</a>
	 *   
	 *   @example InstantPitch
	 */
	public float getFundFreq(AudioSource x){
	  // compute FFT
	  FFT f = new FFT(x.bufferSize(), x.sampleRate());
	  f.window(FFT.HAMMING);
	  f.forward(x.mix);
	  final int ss = f.specSize();  // spectrum size
	  
	  // downsample the spectrum by a factor of 2
	  FFT fd2 = new FFT(x.bufferSize(), x.sampleRate() / 2);
	  fd2.window(FFT.HAMMING);
	  fd2.forward(x.mix);
	  
	  // downsample the spectrum by a factor of 3
	  FFT fd3 = new FFT(x.bufferSize(), x.sampleRate() / 3);
	  fd3.window(FFT.HAMMING);
	  fd3.forward(x.mix);
	  
	  // multiply the trhee spectra together
	  float[] y = new float[ss];
	  for(int i=0; i<ss-1; i++){
	    y[i] = f.getBand(i) * fd2.getBand(i) * fd3.getBand(i);
	  }
	  
	  // find position of the maximum peak in the resulting spectrum
	  float f0Val = -1;
	  int f0Pos = -1;
	  for(int i=0; i<ss; i++){
	    if(y[i] > f0Val){
	      f0Pos = i;
	      f0Val = y[i];
	    }
	  }
	  
	  // convert the index into frequency
	  float f0 = f.indexToFreq(f0Pos);
	  return(f0);
	}//end getFundFreq

	/**
	 * Find peaks within frequency spectrum
	 */
	private void findFreqPeaks(FFT f, ArrayList<Integer> peaks){
	  int ss = f.specSize();
	  
	  float currentDiff = 0;
	  int posOk = -2;
	  final float percRatio = (float) 0.8;
	  
	  float mf = 0;
	  for(int i=0; i<ss; i++){
	    mf += f.getBand(i);
	  }
	  mf /= ss;
	//  println(mf);
	  
	  for(int i=0; i<ss-1; i++){
	    currentDiff = f.getBand(i+1) - f.getBand(i);  // discrete derivative
	    if(currentDiff > 0){  // ascent
	      posOk = i;
	    }
	    else if(currentDiff < 0){  // descent
	      if(posOk == i-1){  // after the peak
	        if((f.getBand(i) > mf * (1 + percRatio)) && (f.getBand(i) > f.getBand(posOk) * (1 - percRatio)) && (mf > 0.01)){
//	        if(f.getBand(posOk) < f.getBand(i) * percRatio){
	          peaks.add(new Integer(i));
	        }
	      }
	    }
	  }
	}//end findFreqPeaks
	
	
	/**
	 * The irregularity of a spectrum is the degree of variation of the successive peaks of the spectrum. It is computed as the sum of the square of the difference in amplitude between adjoining partials.
	 * 
	 *  @param x
	 * 			The time-domain signal to be analyzed
	 * @return irr
	 * 			Spectrum Irregularity coefficient
	 */
	public float getIrregularity(AudioSource x){
	  // compute FFT
	  FFT f = new FFT(x.bufferSize(), x.sampleRate());
	  f.window(FFT.HAMMING);
	  f.forward(x.mix);
	  
	  float irr = 0;  // result
	  float a = 0;
	  float b = 0; 
	  float num = 0;
	  float den = 0;
	  
	  // compute the spectrum peaks
	  ArrayList<Integer> peaks = new ArrayList<Integer>();
	  findFreqPeaks(f, peaks);
	  final int N = peaks.size();
	  
	  // compute the irregularity between adjacent peaks
	  int idx1, idx2;  // peaks indexes
	  for(int i=0; i<N-1; i++){
	    idx1 = (int) peaks.get(i);
	    idx2 = (int) peaks.get(i+1);
	    a = f.getBand(idx1);
	    b = f.getBand(idx2);
	    
	    num += (Math.pow(a - b, 2));
	    den += Math.pow(a, 2);
	  }
	  
	  if(den != 0){
		  irr = num / den;
	  }
	  else{
		  irr = -1;
	  }
	  return(irr);
	}//end getIrregularity
	
	
	/**
	 * The kurtosis gives a measure of the flatness 
	 * of the spectrum around its mean value. 
	 * It is computed from the 4th order of the moment.
	 * 
	 *  @param x
	 * 			The time-domain signal to be analyzed
	 *  @return kurt
	 *  		Spectral Kurtosis coefficient
	 */

	public float getKurtosis(AudioSource x){
	  // compute FFT
	  FFT f = new FFT(x.bufferSize(), x.sampleRate());
	  f.window(FFT.HAMMING);
	  f.forward(x.mix);
	  int ss = f.specSize();  // spectrum size
	  
	  float kurt = 0;  // result
	  float xi = 0;
	  float fi = 0;
	  
	  // compute centroid and spread
	  float sc = getCentroid(x);  // spectrum centroid
	  float sigma = (float) Math.sqrt(getSpread(x));  // spectrum standard deviation
	  
	  if(sigma != 0){
		  // compute kurtosis
		  for(int i=0; i<ss-1; i++){  // for each frequency band...
			  xi = f.getBand(i);
			  fi = f.indexToFreq(i);
			  kurt += (xi * (Math.pow((fi - sc), 4)));  // Why 4? By definition...
		  }

		  kurt /= Math.pow(sigma, 4);
	  }
	  else{
		  kurt = -1;
	  }
	  return(kurt);
	}//end getKurtosis
	
	
	/**
	 * Mel-Frequency Cepstral Coefficients (MFCC) originate 
	 * from automatic speech recognition.
	 * The MFCCs are to some extent created according to 
	 * the principles of the human auditory system, but also 
	 * to be a compact representation of the amplitude 
	 * spectrum and with considerations of the computational 
	 * complexity.
	 * 
	 * @param x
	 * 			The time-domain signal to be analyzed
	 * @param index
	 * 			The order of the MFCC to be computed (usually 1 ² <code>index</code> ² 13)
	 * @return melCoeff
	 * 			The <code>index</code>-th order Mel Cepstral coefficient
	 * @see <a href="http://www.speech-recognition.de/matlab-files.html">Mel filterbank algorithm</a>
	 */

	public float getMFCC(AudioSource x, int index){
	  // compute FFT
	  FFT f = new FFT(x.bufferSize(), x.sampleRate());
	  f.window(FFT.HAMMING);
	  f.forward(x.mix);
	  final int N = f.specSize();  // spectrum size
	  
	  float fs = x.sampleRate();  // sampling frequency
	  final int nofChannels = 22;  // number of mel channels
	  
	  // compute resolution etc
	  final float df = fs / N;  //  frequency increment on linear scale
	  final float Nmax = N; //  maximum fft index
	  final float fmax = fs/2; //  maximum frequency
	  final float melmax = freq2mel(fmax); //  maximum mel frequency
	  final float melinc = melmax / (nofChannels + 1); //  frequency increment on mel scale
	  
	  //  compute center frequencies on mel scale
	  float[] melcenters = new float[nofChannels];
	  for(int i=0; i<nofChannels; i++){
	    melcenters[i] = (i + 1) * melinc;
	  }
	  
	  // compute center frequencies in linear scale [Hz]
	  float[] fcenters = new float[nofChannels];
	  for(int i=0; i<nofChannels; i++){
	    fcenters[i] = mel2freq(melcenters[i]);
	  }
	  
	  // compute indexes of center frequencies
	  float[] indexcenter = new float[nofChannels];
	  for(int i=0; i<nofChannels; i++){
	    indexcenter[i] = Math.round(fcenters[i] / df);
	  }
	  
	  float[] indexstart = new float[nofChannels];  // compute start indices of windows
	  float[] indexstop = new float[nofChannels];  // compute stop indices of windows
	  float[] idxbw = new float[nofChannels];  // compute bandwidth (number of indices per window)
	  for(int i=0; i<nofChannels; i++){
	    if(i == 0){
	      indexstart[i] = 0;
	    }
	    else if(i == nofChannels - 1){
	      indexstop[i] = Nmax-1;
	    }
	    else{
	      indexstart[i] = indexcenter[i-1];
	      indexstop[i] = indexcenter[i+1];
	    }
	    idxbw[i] = (indexstop[i] - indexstart[i]) + 1;
	  }//end for
	  
	  // signal spectrum in mel components
	  float[] melSpectrum = new float[N];
	  for(int i=0; i<N; i++){
	    melSpectrum[i] = 0;
	  }
	  for(int i=0; i<N; i++){
	    int currMelChannel = 0;
	    for(int j=0; j<nofChannels-1; j++){
	      if(i >= indexstart[j] && i < indexstart[j+1]){
	        currMelChannel = j;  // current mel channel
	      }
	    }
	    
	    float Si = f.getBand(i);
	    float Nmel = idxbw[currMelChannel];
	    float currMelIndex = indexstop[currMelChannel] - i;
	    melSpectrum[i] += Si * (2.0 / (Nmel - 1.0) * (((Nmel - 1.0) / 2.0) - Math.abs(currMelIndex - ((Nmel - 1.0) / 2.0))));
	  }//end for
	  
	  // take the log of the mel frequency components
	  for(int i=0; i<N; i++){
	    if(melSpectrum[i] > 0){  // avoid -Infinity
	      float temp = (float) Math.log(melSpectrum[i]);
	      melSpectrum[i] = temp;
	    }
	    else{
	      melSpectrum[i] = 0;
	    }
	  }
	  
	  // perform DCT on the log values
	  float melCoeff = 0;
	  for(int i=0; i<N-1; i++){
	    melCoeff += melSpectrum[i] * Math.cos((Math.PI/N) * (i + 0.5) * index);
	  }
	  
	  return(melCoeff);
	}//end getMFCC
	
	
	/**
	 * The spectrum spread is also called instantaneous 
	 * bandwidth; it is defined as the second central 
	 * moment (i.e. the variance) of spectrum.
	 * 
	 *  @param x
	 * 			The time-domain signal to be analyzed
	 *  @return ass
	 *  		Audio Spectrum Spread coefficient
	 */
	public float getSpread(AudioSource x){
	  // compute FFT
	  FFT f = new FFT(x.bufferSize(), x.sampleRate());
	  f.window(FFT.HAMMING);
	  f.forward(x.mix);
	  final int Nft = f.specSize();  // spectrum size
	  
	  float deltaF = x.sampleRate() / Nft;  // frequency interval between two FFT bins
	  int kLow = (int) Math.floor(62.5 / deltaF);  // preventing disproportionate weight below this frequency
	  
	  float asc = getCentroid(x);  // compute spectrum centroid
	  
	  float ass = 0;  // result
	  float num = 0;
	  float den = 0;
	  
	  // compute first component
	  float Pk = 0;
	  if(kLow > 0){
	    for(int i=0; i<kLow; i++){  // for each band
	      Pk += f.getBand(i);
	    }
	    num += (31.25 - asc) * Pk;  // 31.25 Hz = fLow
	    den += Pk;
	  }
	  
	  // compute the remaining components
	  float fk = 0;
	  for(int i=kLow; i<Nft/2; i++){
	    Pk = f.getBand(i);
	    fk = f.indexToFreq(i);
	    num += (fk - asc) * Pk;
	    den += Pk;
	  }
	  
	  ass = num / den;
	  return(ass);
	}//end getSpread
	
	
	/**
	 * Short Time Energy (STE) expresses the instant energy of each frame, and is commonly known as volume.
	 * 
	 * A speech audio sample will give large variations in
	 * the STE-values. Music has much shorter pauses, or no 
	 * pauses, and will therefore have a more constant 
	 * STE-level.
	 *   @param x
	 * 			The time-domain signal to be analyzed
	 *   @return ste
	 *   		Short Time Energy coefficient
	 */
	public float getSTE(AudioSource x){
	  final int L = x.bufferSize();
	  int ste = 1;  // result
	  float a;
	  
	  // compute STE
	  for(int i=0; i<L; i++){
	    a = x.mix.get(i);
	    ste += (a * a);
	  }
	  
	  return(ste);
	  
	}//end getSTE
	
	
	
	/**
	 * The zero-crossing rate is the rate of sign-changes
	 * along a signal, i.e., the rate at which the signal 
	 * changes from positive to negative or back.
	 * 
	 *  @param x
	 * 			The time-domain signal to be analyzed
	 * @return z
	 * 			Amount of sign changes within the input signal
	 */
	public int getZCR(AudioSource x){
	  final int L = x.bufferSize();
	  int z = 0;  // result
	  float a, b;
	  
	  // compute ZCR
	  for(int i=0; i<L-1; i++){
	    a = x.mix.get(i);
	    b = x.mix.get(i+1);
	    if(a*b < 0){  // different signs
	      z++;
	    }
	  }
	  
	  return(z);
	}//end getZCR
	
	
	/**
	 * The skewness gives a measure of the asymmetry of the 
	 * spectrum around its mean value. It is computed from 
	 * the 3rd order of the moment.
	 *  
	 *  @param x
	 * 			The time-domain signal to be analyzed
	 *  @return skew
	 *  		Spectrum Skewness
	 */
	public float getSkewness(AudioSource x){
	  // compute FFT
	  FFT f = new FFT(x.bufferSize(), x.sampleRate());
	  f.window(FFT.HAMMING);
	  f.forward(x.mix);
	  int ss = f.specSize();  // spectrum size
	  
	  float skew = 0;  // result
	  float xi = 0;
	  float fi = 0;
	  
	  // compute spectrum centroid and "standard deviation"
	  float sc = getCentroid(x);  // spectrum centroid
	  float sigma = (float) Math.sqrt(getSpread(x));  // spectrum standard deviation
	  
	  // compute skewness
	  for(int i=0; i<ss-1; i++){  // for each frequency band...
	    xi = f.getBand(i);
	    fi = f.indexToFreq(i);
	    skew += (xi * (Math.pow((fi - sc), 3)));  // Why 3? By definition...
	  }
	  
	  skew /= Math.pow(sigma, 3);
	  return(skew);
	}//end getSkewness
	
	
	/**
	 * The roughness of a signal is related to the 
	 * beating phenomena between the peaks of its 
	 * spectrum.
	 * 
	 * An estimation of the total roughness can be 
	 * obtained by computing the peaks of the 
	 * spectrum, and taking the average of all the 
	 * dissonance between all possible pairs of 
	 * peaks.
	 * 
	 *  @param x
	 * 			The time-domain signal to be analyzed
	 *  @return r
	 *  		Roughness coefficient
	 *  
	 *  @example CharlestonDetection
	 */

	public float getRoughness(AudioSource x){
	  // compute FFT
	  FFT f = new FFT(x.bufferSize(), x.sampleRate());
	  f.window(FFT.HAMMING);
	  f.forward(x.mix);
	  ArrayList<Integer> peaks = new ArrayList<Integer>();  // spectrum peaks
	  
	  float a = 0;
	  float b = 0;
	  int aPos, bPos;
	  float Amin, Amax, fmin, fmax, X, Y, Z;
	  
	  // constant coefficients
	  final float b1 = (float) 3.5;
	  final float b2 = (float) 5.75;
	  final float s1 = (float) 0.0207;
	  final float s2 = (float) 18.96;
	  
	  float r = 0;  // result
	  
	  // find frequency peaks
	  findFreqPeaks(f, peaks);
	  final int K = peaks.size();
	  
	  for(int i=0; i<K-1; i++){
	    for(int j=i+1; j<K; j++){
	      aPos = (Integer) peaks.get(i);
	      bPos = (Integer) peaks.get(j);
	      a = f.getBand(aPos);
	      b = f.getBand(bPos);
	      Amin = Math.min(a, b);  // minimum frequency component
	      Amax = Math.max(a, b);  // maximum frequency component
	      fmin = f.indexToFreq(aPos);
	      fmax = f.indexToFreq(bPos);
	      
	      final float s = (float) (0.24 / (s1 * fmin + s2));
	      
	      X = Amin * Amax;
	      Y = (2 * Amin) / (Amin + Amax);
	      Z = (float) (Math.exp(-b1 * s * (fmax - fmin)) - Math.exp(-b2 * s * (fmax - fmin)));
	      
	      r += (Math.pow(X, 0.1) * 0.5 * (Math.pow(Y, 3.11)) * Z);
	    }//end for
	  }//end for
	  
	  r /= K;
	  return(r);
	}//end getRoughness
	
	
}//end class


