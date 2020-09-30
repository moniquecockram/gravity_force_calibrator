// NCal control code
// PHYS6701 S2 2020
// Monique Cockram and Ryan Husband


//int pins for the motor and photogate
//int motor =  
int sensor = 19;
//bool sensorValue; 


void setup() {
  // put your setup code here, to run once:

  Serial.begin(9600);

//initialise pins
  //pinMode(motor, OUTPUT)

  pinMode(sensor, INPUT);

 
//do a slow speed up of the wheel
//write(motor, )
//delay(1000)
//write(motor, slightly faster)
//delay(1000)
//etc


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



  //sensorValue = digitalRead(sensor);
 
  Serial.println(digitalRead(sensor));
 
  delay(500);





}


//figure out how to do a shut down to slowly reduce the speed?
