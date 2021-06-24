E = 210e9; q = 0; P = 2000; M = 0; density = 7850;b = 0.05; h = 0.01; I =b*h*h*h/12;

le = 1;
Element_no = 20;
elem_len = le/Element_no;
Nodes = Element_no + 1;
DOF = Nodes*2;
f = zeros(4,1,Element_no);
k = zeros(4,4,Element_no);
Mass = zeros(4,4,Element_no);
%Element stiffness matrix and load vector
for element = 1:1:Element_no
    k(:,:,element) = Beam_stiff_matrix(E,I,[(element-1)*elem_len,element*elem_len]);
    Mass(:,:,element) = Beam_Mass_matrix(b*h,density,elem_len);
    f(:,:,element) = Load_vec_matrix(q,[(element-1)*elem_len,element*elem_len]);
end

K_final = zeros(DOF,DOF);
M_final = zeros(DOF,DOF);
K_final([1:4],[1:4]) = k(:,:,1);
M_final([1:4],[1:4]) = Mass(:,:,1);
F_final = zeros(DOF,1);
F_final(1:4) = f(:,:,1);
%Assembling
for element = 1:1:(Element_no-1)
    K_final([(element*2)+1:(element*2)+4],[(element*2)+1:(element*2)+4]) = K_final([(element*2)+1:(element*2)+4],[(element*2)+1:(element*2)+4]) + k(:,:,element+1);
    M_final([(element*2)+1:(element*2)+4],[(element*2)+1:(element*2)+4]) = M_final([(element*2)+1:(element*2)+4],[(element*2)+1:(element*2)+4]) + Mass(:,:,element+1);
    F_final([(element*2)+1:(element*2)+4]) = F_final([(element*2)+1:(element*2)+4]) + f(:,:,element+1);
end

% k1 = Beam_stiff_matrix(E,I,[0,le]);
% f = Load_vec_matrix(q,[0,le]);
F_final(DOF-1) = F_final(DOF-1) + P;
F_final(DOF) = F_final(DOF) + M;
% f(3) = f(3) + P;
% f(4) = f(4) + M;

%Solving equations for deflection ( Essential Boundary conditions : w1 = 0
%and (dw/dx)1 = 0)
W = zeros(DOF,1);

Kreduced = K_final([3:DOF],[3:DOF]);
Freduced = F_final([3:DOF]);
K_inv = inv(Kreduced);
wreduced = K_inv*Freduced;

W(3:DOF) = wreduced;
Reaction = K_final*W;
disp("SUCCESS");

%Natural Frequency
freq_sqr = eig(K_final,M_final);
freq = (1/(2*pi))*sqrt(freq_sqr);
%Mass_mat = Beam_Mass_matrix(1,density,le);


%Plotting the solution
 w_final = zeros(Element_no*length([-1,0.1,1]),1);
 w_actual = zeros(Element_no*length([-1,0.1,1]),1);
 X_vals = zeros(Element_no*length([-1,0.1,1]),1);
 i = 1;
 for element = 1:1:Element_no
 for e = -1 : 0.1 : 1
     N1 = 0.25*(1-e)*(1-e)*(2+e);
     N2 = 0.25*(1-e)*(1-e)*(1+e);
     N3 = 0.25*(1+e)*(1+e)*(2-e);
     N4 = 0.25*(1+e)*(1+e)*(e-1);
     w_final(i) = N1*W(2*element - 1) + elem_len*0.5*N2*W(2*element) + N3*W(2*element + 1) + elem_len*0.5*N4*W(2*element +2);
     X_fin = (1+e)*0.5*element*elem_len + (1-e)*0.5*(element-1)*elem_len;
     X_vals(i) = X_fin;
     w_actual(i) = (P*X_fin*X_fin/(6*E*I))*(3*le - X_fin);
     i = i+1;
 end
 end

plot(X_vals,-1*w_final,'*',X_vals,-1*w_actual);
legend('FEM results','Analytical Solution');
title('20 Element Solution');
xlabel('X(m)');ylabel('Deflection(m)');
% plot(linspace(1,DOF,42),freq(1:DOF),'*',linspace(1,DOF,42),freq(1:DOF));
% xlabel('DOF');ylabel('Natural Frequency(Hz)');
% x = linspace(0,10,21);
% plot(x,-1*w_final,'blue');
% plot(x,-1*w_actual,'red');
% title('Analytical Solution'); xlabel('x(m)');ylabel('Deflection(m)');



