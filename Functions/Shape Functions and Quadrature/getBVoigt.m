function [Bvoigt] = getBVoigt(Mesh, B)
% return the B matrix in Voigt form
% Input: NB matrix
% adapted from GitHub GCMLab-FEM, getNv.m
% Acknowledgements: Matin Parchei Esfahani

nsd = Mesh.nsd;

switch nsd
    case 1
        Bvoigt = B;
    case 2
        if Mesh.field == 'p'
            Bvoigt = B;
        else
            Bvoigt = zeros(3, Mesh.nsd*Mesh.nne);

            Bvoigt(1,1:2:end) = B(1,:);
            Bvoigt(2,2:2:end) = B(2,:);
            Bvoigt(3,1:2:end) = B(2,:);
            Bvoigt(3,2:2:end) = B(1,:);
        end
end

end