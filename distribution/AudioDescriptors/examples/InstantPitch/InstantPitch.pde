/**
 * Instant Pitch
 * by Davide Totaro.
 *  
 * This sketch demonstrates how to use the <code>getFundFrequency</code> method 
 * of <code>AudioDescriptor</code>, a class totally compatible with standard 
 * Minim library.
 * The <code>getFundFrequency</code> method allows you to extract the Fundamental 
 * Frequency from the <code>mix/code> channel of an<code>AudioPlayer</code> 
 * object. The returned value is measured in Hertz [Hz].
 * 
 * Every AudioDescriptor function receives as input an <code>AudioSource</code> 
 * object, like <code>AudioPlayer</code>.
 */

import ddf.minim.*;
import ddf.minim.analysis.*;

import audio.descriptors.*;

// define Minim objects
Minim minim;
AudioInput in;
FFT f;
int N;  // number of spectrum components

// define AudioDescriptor object
AudioDescriptor ad;

// declare and initialize constant values
final int mulFactor = 11;

void setup(){
  size(512, 200, P2D);

  minim = new Minim(this);
  minim.debugOn();
  
  // get a line in from Minim, default bit depth is 16
  in = minim.getLineIn(Minim.STEREO, 512);
  
  // create and set the FFT object
  f = new FFT(in.bufferSize(), in.sampleRate());
  f.window(FFT.HAMMING);
  N = f.specSize();
  
  // create Audio Descriptor through constructor
  ad = new AudioDescriptor();
  
}//end setup

void draw(){
  background(0);
  stroke(255);
  
  // compute the frequency components of the input signal
  f.forward(in.mix);
  
  // draw the spectrum
  for(int i = 0; i < N - 1; i++){
    line(i*2, height - 5, i*2, height - 5 - f.getBand(i) * mulFactor);
  }
  
  // compute the Fundamental Frequency of the input mix channel.
  // this value is measured in Hertz [Hz].
  // show it on the spectrum with a different color.
  float f0 = ad.getFundFreq(in);
  stroke(11,11,222);
  int k0 = f.freqToIndex(f0);  // convert the frequency from Hertz into FFT index
  line(k0*2, height - 11, k0*2, height - 11 - f.getBand(k0) * mulFactor);
}//end draw


void stop(){
  in.close();
  minim.stop();
  
  super.stop();
}
