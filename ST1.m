function Result = ST1(R,G,B)
 Result= G-R>0;

% Alternatively, other indices such as the formula below can be used for ST1, depended on the experimental conditions: 

% Result= G>R&G>B;
% Result= (2*G-R-B)>0;

end
