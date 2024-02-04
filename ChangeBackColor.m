function Res = ChangeBackColor(I,r,g,b)
    R = I(:,:,1);
    G = I(:,:,2);
    B = I(:,:,3);
    Ind_R = find(R==0);
    R(Ind_R) = r;
    Ind_G = find(G==0);
    G(Ind_G) = g;
    Ind_B = find(B==0);
    B(Ind_B) = b;
    Res(:,:,1) = R;
    Res(:,:,2) = G;
    Res(:,:,3) = B;
end