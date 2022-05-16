classdef DeltaSigmaModulator < hgsetget
    
    properties
        Sigma          = 0;
        Threshold      = 1;
        PreviousOutput = 0;
        Oversampling   = 1;
    end
    
    methods
    
         % Constructor
        function Object = DeltaSigmaModulator(varargin)
            
            for v = 1:2:length(varargin)
                Property = varargin{v};
                Value = varargin{v+1}; 
                set(Object,Property,Value);                
            end
            
        end
        
        % Function 'set'
        function Object = set(Object,varargin)
            
            Properties = varargin(1:2:end);
            Values = varargin(2:2:end);            
            for n = 1:length(Properties)                
                [is, Property] = isproperty(Object,Properties{n}); 
                if is      
                    switch Property
                        case 'Oversampling'
                            if Values{n} ~= floor(Values{n})
                                error('Oversampling factor must be an integer value.')
                            end
                        case 'PreviousOutput'
                            if Values{n} < 0 || Values{n} > 1
                                error('Previous output must be a scalar between 0 and 1.')
                            end
                    end                    
                    Object.(Property) = Values{n};       
                else
                    error('Property "%s" not supported !',Properties{n});
                end                
            end
            
        end
        
        % Function 'get'
        function Value = get(varargin)
            
            switch nargin                
                case 1                    
                    Value = varargin{1};                    
                otherwise                    
                    Object = varargin{1};
                    [is, Property] = isproperty(Object,varargin{2});
                    if is                        
                        Value = Object.(Property);
                    else
                        error('Property "%s" not supported !',varargin{2});
                    end
                    
            end
            
        end
        
        % Function 'isproperty'
        function [is, Property] = isproperty(Object,Property)
            
            Properties = fieldnames(Object); 
            [is, b] = ismember(lower(Property),lower(Properties));
            Property = Properties{b};
            
        end
   
        % Function 'update'
        function varargout = update(Object,Input)
                        
            I = size(Input);
            
            % Oversampling by holding
            if I(1) >= I(2)                
                Input = repmat(Input',Object.Oversampling,1);
                Input = reshape(Input,[],1);
                S1 = I(1) * Object.Oversampling;
                S2 = 1;
            else                
                Input = repmat(Input,Object.Oversampling,1);
                Input = reshape(Input,1,[]);
                S1 = 1;
                S2 = I(2) * Object.Oversampling;
            end
            
            % Modulation
            Output = zeros(S1,S2);            
            for s = 1:max(S1,S2)
                Delta =  Input(s) - Object.PreviousOutput;
                Object.Sigma = Object.Sigma + Delta;
                Output(s) = ge(Object.Sigma,Object.Threshold);
                Object.PreviousOutput = Output(s);
            end
            
            % Outputs definition
            switch nargout
                case 1
                    varargout{1} = Output;
                case 2
                    varargout{1} = Input;
                    varargout{2} = Output;
            end
            
        end
        
    end
        
end



