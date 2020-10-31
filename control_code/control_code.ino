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


int BPin = 6; //backwards motor in1
int FPin = 7; //forwards motor in2

int speedPin = 19; //enA 



float chosen_period = 15.0; //time in seconds for one rotation


float period; //the period of rotation of the wheel
float num_holes = 5; //the number of holes around the circumference of the wheel
unsigned long prev_time = 0;
unsigned long current_time = 0;
float time_taken = 0;
int speed_val;


void setup() {
  // put your setup code here, to run once:

  Serial.begin(9600);
  

//initialise pins

  pinMode(13, OUTPUT);
  pinMode(gatePin, INPUT);


  pinMode(BPin, OUTPUT);
  pinMode(FPin, OUTPUT);
  pinMode(speedPin, OUTPUT);
  

  attachInterrupt(digitalPinToInterrupt(gatePin), timing_function, RISING);


  digitalWrite(BPin,LOW);
  digitalWrite(FPin,HIGH);
  speed_val = 150; //an initial speed value
  analogWrite(speedPin, speed_val); //previous projects showed that this is lowest speed to overcome friction


}


void timing_function() {
  //interupt function to execute when the photogate changes to HIGH
  current_time = millis();
  time_taken = current_time - prev_time;
  prev_time = current_time;
  if ((period > chosen_period) && (speed_val < 255)) {
    speed_val = speed_val + 1; //increasing speed_val if too slow
    analogWrite(speedPin, speed_val);
    }
    
  if ((period < chosen_period) && (speed_val > 128)) { //128 is the lowest speed before its too slow for friction
    speed_val = speed_val - 1;
    analogWrite(speedPin, speed_val);
    }
  }




void loop() {
  // put your main code here, to run repeatedly:
 
   //Read photogate data
  val = digitalRead(gatePin);
  //Print photogate data in binary
  //Serial.print(val); //printing photogate data
  //Serial.print("\n");

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
  Serial.print("___");
  period = time_taken*num_holes/1000; //in seconds
  Serial.print(period);
  Serial.print("___");
  delay(10); 
  Serial.print(speed_val);
  Serial.print("\n");
  
}
