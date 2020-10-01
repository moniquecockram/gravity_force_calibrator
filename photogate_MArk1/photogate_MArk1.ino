/*<Code stolen from <https://github.com/sparkfun/Photo_Interrupter_Breakout> used on https://www.youtube.com/watch?v=jNMb2Az-hIY&ab_channel=SparkFunElectronics
  Simple Sketch for getting started with Photo Interrupter Breakout */

//Declare libraries
#include <SoftwareSerial.h>

//Set pins
const int gatePin = 10;
int val;

void setup()
{
  //Set pin modes
  pinMode(13, OUTPUT);
  pinMode(gatePin, INPUT);
  Serial.begin(9600);
}

void loop()
{
  //Read photogate data
  val = digitalRead(gatePin);
  //Print photogate data in binary
  Serial.print(val);
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

delay(10);
}
