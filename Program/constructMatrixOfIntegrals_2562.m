function Bg=constructMatrixOfIntegrals_2562(GradientOrientations, order, delta)

g_=GradientOrientations;
%UnitVectors;
%g2=g([1:321],:);
g2=BuildSphere(4); % 2562 directions
%g = sphere_tri('ico',4,1);
%g2=g.vertices;

%load sphere_triangulation_2562
%g2=cen; % centroids

for k=1:size(g_,1)
    c=1;
    for i=0:order
        for j=0:order-i
            Bg(k,c)=0;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %integr=1:size(g2,1);
            %Bg(k,c)=exp(-delta*(g_(k,:)*g2').^2)*((g2(:,1).^i).*(g2(:,2).^j).*(g2(:,3).^(order-i-j)).*areas');
            Bg(k,c)=exp(-delta*(g_(k,:)*g2').^2)*((g2(:,1).^i).*(g2(:,2).^j).*(g2(:,3).^(order-i-j)));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Bg(k,c)=4*pi*Bg(k,c)/size(g2,1);
            %Bg(k,c)=Bg(k,c).*areas'; % multiply by the areas of the spherical triangles
            
            c=c+1;
        end
    end
end
%fprintf(1,'Done\n');