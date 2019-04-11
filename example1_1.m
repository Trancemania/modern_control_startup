s = tf('s');
tf = (2*s+4)/(s^2*(s^2+3*s+6));
num = cell2mat(tf.numerator(1,1));
den = cell2mat(tf.denominator(1,1));