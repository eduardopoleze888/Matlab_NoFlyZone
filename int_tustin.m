function yn = int_tustin(un, u_n1, yn_1, Td)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    yn = (Td/2)*(un + u_n1) + yn_1;

end

