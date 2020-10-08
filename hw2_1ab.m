A = [4 -1; 1 0];

% plot ||b||=1, ||Ab||=1, 2-norm
% no need to generalize here, we know sol is x^2+y^2=1
% and for Ab=1 is (4x-y)^2+x^2=1

% make implicit plot (=0)
fimplicit(@(x,y) x.^2 + y.^2 - 1);
hold on
fimplicit(@(x,y) (4*x-y).^2 + x.^2 - 1);
legend(["||b||_2=1" "||Ab||_2=1"], 'Location', 'northwest')

% make the other implicit plot
% infinity-norm
%hold off;
fimplicit(@(x,y) max(abs(x),abs(y)) - 1);
hold on
fimplicit(@(x,y) max(abs(4*x-y), abs(x)) - 1);
legend(["||b||_\infty=1" "||Ab||_\infty=1"], 'Location', 'northwest')

% make implicit plot (=0) for 4-norm
fimplicit(@(x,y) (x.^4 + y.^4).^(1/4) - 1);
hold on
fimplicit(@(x,y) ((4*x-y).^4 + x.^4).^(1/4) - 1);
legend(["||b||_4=1" "||Ab||_4=1"], 'Location', 'northwest')