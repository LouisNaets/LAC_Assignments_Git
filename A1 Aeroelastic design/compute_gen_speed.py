v_rated = 11.2
R = 92.4
gear_ratio = 50
rotor_speed = 9.6

omega = rotor_speed * 2 * 3.14159 / 60          #radians /s
TSR = omega*R/ v_rated

genspeed = TSR * gear_ratio         #411

print(f'TSR = {TSR}')
print(f'genspeed = {genspeed}')