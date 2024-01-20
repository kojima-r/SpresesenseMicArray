import serial

readSer = serial.Serial('/dev/ttyUSB0',115200, timeout=3)
#c = readSer.read() # 1 byte
#string = readSer.read(10) # 10 byte
while True:
    line = readSer.readline() # 1 line (upto '\n')
    print("Read Serial:")
    print(line)

readSer.close()
