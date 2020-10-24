// NCal control code
// PHYS6701 S2 2020
// Monique Cockram and Ryan Husband


/*<Code stolen from <https://github.com/sparkfun/Photo_Interrupter_Breakout> 
 * used on https://www.youtube.com/watch?v=jNMb2Az-hIY&ab_channel=SparkFunElectronics
  Simple Sketch for getting started with Photo Interrupter Breakout */

//Declare libraries
#include <SoftwareSerial.h>

//set pins for the motor and photogate
const int gatePin = 10;
int val;



unsigned long prev_time = 0;
unsigned long current_time = 0;
unsigned long time_taken = 0;



void setup() {
  // put your setup code here, to run once:

  Serial.begin(9600);
  

//initialise pins

  pinMode(13, OUTPUT);
  pinMode(gatePin, INPUT);

 


  attachInterrupt(digitalPinToInterrupt(gatePin), timing_function, RISING);
 
//do a slow speed up of the wheel
//write(motor, )
//delay(1000)
//write(motor, slightly faster)
//delay(1000)
//etc


}


void timing_function() {
  //interupt function to execute when the photogate changes to HIGH
  current_time = millis();
  time_taken = current_time - prev_time;
  prev_time = current_time;
  //Serial.print("interupted");
  }




void loop() {
  // put your main code here, to run repeatedly:


// set wheel speed to chosen value
// read light gate output data (should be a time between last peak?
// calculate rotation speed speed = (C/number of holes)/time = (2*pi*r/50)/time
// adjust rotation speed accordingly (if loop?)
      // if speed < chosen value
            // set_speed += some small value (related to the error from chosen speed? e.g. chosen speed - measurered speed)
      // if speed > chosen value
            // set_speed -= some small value
     // else do nothing


 
   //Read photogate data
  val = digitalRead(gatePin);
  //Print photogate data in binary
  Serial.print(val); //printing photogate data
  Serial.print("\n");

  //If photogate detects something, turn Teensy LED on
  if (val == HIGH) {
    digitalWrite(LED_BUILTIN, HIGH);
  }

  //Leave LED off otherwise
  else
  {
    digitalWrite(LED_BUILTIN, LOW);
  }

  Serial.print(time_taken);
  Serial.print("\n");
  delay(10); //make sure this delay is short enough that it wont miss a hole, or even just remove it




}


//figure out how to do a shut down to slowly reduce the speed?
