function [Model] = get_model(Nx, model_choice)

switch model_choice
    case 'BS' % Black Scholes Model
%         Model = Black_Scholes_Model(Nx);
        Model = BS_Model1();
    case 'md_RBF' % Mahalanobis distance RBF
        Model = RBF1_Model(Nx);
        
    case 'gaus_RBF' % Gaussian RBF
        Model = get_RBF2(Nx);
end

end