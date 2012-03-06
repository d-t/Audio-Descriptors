/**
 * Charleston Detection
 * by Davide Totaro.
 *  
 * This sketch demonstrates how to use the <code>getRoughness</code> method 
 * of <code>AudioDescriptor</code>, a class totally compatible with standard 
 * Minim library.
 * The <code>getRoughness</code> method allows you to extract the Roughness 
 * coefficient from the <code>mix/code> channel of an<code>AudioPlayer</code> 
 * object. The returned value is contained between 0 and 1.
 * 
 * Every AudioDescriptor function receives as input an <code>AudioSource</code> 
 * object, like <code>AudioPlayer</code>.
 */

import ddf.minim.*;
import audio.descriptors.*;

// define Minim objects
AudioPlayer player;
Minim minim;

// define AudioDescriptor object
AudioDescriptor ad;

// declare and initialize constant values
final float mulFactor = 40;
final int rad = 111;  // radius
final float thresholdValue = 0.055;

void setup(){
  size(512, 200, P2D);

  minim = new Minim(this);
  // load a file, give the AudioPlayer buffers that are 1024 samples long
  player = minim.loadFile("groove.mp3");
  
  // play the file
  player.play();
  
  // create Audio Descriptor through constructor
  ad = new AudioDescriptor();
  
  // drawing parameters
  noStroke();
  fill(222, 111, 0);
  smooth();
}//end setup

void draw(){
  background(0);
  
  // compute the Roughness coefficient of the player mix channel and
  // detect the charleston sound.
  // since this value will be between 0 and 1, we use an exponential
  // function in order to better show its effects
  float rough = ad.getRoughness(player);
  if(rough > thresholdValue){  // Roughness coefficient greater than predefined threshold
    ellipse(width/2, height/2, rad + exp(rough*mulFactor), rad + exp(rough*mulFactor));
  }
  else{
    ellipse(width/2, height/2, rad, rad);
  }
}//end draw

void stop(){
  player.close();
  minim.stop();
  
  super.stop();
}
