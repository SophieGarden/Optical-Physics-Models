/* 
 functions:
 1 move the motor according to the wavelength in a set speed 
 2 read from PMT/Photodiode->Boxcar->Auduino for the fluorescence and pump intensity
 
 */

#include <Stepper.h>
//#include <stdlib.h>

// Stuff for serialEvent()
String inputString = "";
boolean stringComplete = false;
int inputsteps = 0;
float timeD = 0;
//char test[] = "";

const int stepsPerRevolution = 200;  // the number of steps per revolution for your motor
long STEPS = 0;  //steps>0 clockwise; steps<0, anticlockwise  the number of steps to turn the motor - positive to turn one direction, negative to turn the other (int)
int stepCount = 0;

// initialize the stepper library on pins 8 through 11:
Stepper myStepper(stepsPerRevolution, 8,9,10,11);    

// Analog read pins
const int LIF_pin = A0;    // laser induced fluorescence (PMT-->Boxcar-->HERE)
const int Laser_pin = A3;  // laser intensity (Photodiode-->Boxcar-->HERE)
const int Busy_pin = A2;   // busy signal ('BUSY' from Boxcar)

// Analog read values
int LIF = 0;        // this will hold the analog read value for the LIF
int Laser = 0;      // ditto for the laser intensity

// other
int BusyThreshold = 204;



void setup() {
  Serial.begin(9600);        // open the serial port at 9600 baud
  inputString.reserve(20);   // reserve 200 bytes for inputString
  myStepper.setSpeed(60);    // rpm
}
  

void loop() {
  // decode the input
if (stringComplete) {
 
//******************************************Read Data part*******************************************
    if (inputString == "PMT"){
      
      // wait for busy
      while(analogRead(Busy_pin) < BusyThreshold){
      }
      
      // wait for not busy
      while(analogRead(Busy_pin) > BusyThreshold){
      }  
      
      delayMicroseconds(500);    // wait for things to settle after transitioning from busy to not-busy
      
      // read the LIF and Laser intensities from the boxcars
      LIF = analogRead(LIF_pin);
      Laser = analogRead(Laser_pin);
      
      // send these values to the serial port
      Serial.print(LIF);
      Serial.print(',');
      Serial.print(Laser);
      Serial.println(":");
    }
    
    
//******************************************Motor Part***********************************************
  else if (inputString != "PMT") {

  //Serial.println(inputString);

  char test_as_char[inputString.length()]; // or:  char *test_as_char; .... not sure if you have to initialize the buffer
  inputString.toCharArray(test_as_char, inputString.length()+1); // aotumatically add an null hence length has to add 1
  STEPS = atol(test_as_char);
 // Serial.println(STEPS);
 // timeD=STEPS/200/30*60*1000+500; // wait for motor
 // Serial.println(timeD/1000);
  myStepper.step(STEPS); 
  
 // Serial.flush();

 //delay(timeD); //unit in millisecond

}
  
//****************************************************************************************************
    // clear the string:
    inputString = "";
    stringComplete = false;
  }
}

void serialEvent() {
  //read string for MTR or PMT
  while (Serial.available()) {
    // get the new byte:
    char inChar = (char)Serial.read(); 
    // add it to the inputString:
    // if the incoming character is a newline, set a flag
    // so the main loop can do something about it:
    if (inChar == '\n') {
      stringComplete = true;
       } 
    else{
      inputString += inChar;
    }
  }
 }
