classdef Email
    %EMAIL
    %   Utility to send emails.
    
    properties (~Dependent, Hidden, GetAccess = private)
        address  = ''                       % Email address.
        password = ''                       % Email password.
        server   = ''                       % Email server.
    end
    
    properties (~Dependent)
        recipients = ''                     % Recipients to send email.
        subject    = ''                     % Subject of the email.
        message    = ''                     % Message to be sent.
    end
    
    methods
        function self = Email(varargin)
            %
            % Email instance.
            %
            % INPUTS:
            %   1. address  - email address.
            %   2. password - email password.
            %   3. server   - email server.
            %   4. recipients - recipients to send email.
            %   5. subject  - email subject.
            %   6. message  - email message.
            %
            narginchk(0,6);
            
            if ~isempty(varargin)
                if nargin >= 6
                    self.message = varargin{6};
                end
                
                if nargin >= 5
                    self.subject = varargin{5};
                end
                
                if nargin >= 4
                    self.recipients = varargin{4};
                end
                
                if nargin >= 3
                    self.server = varargin{3};
                end
                
                if nargin >= 2
                    self.password = varargin{2};
                end
                
                if nargin >= 1
                    self.address = varargin{1};
                end
            end
        end
        
        function self = set.address(self, newAddress)
            self.address = newAddress;
        end
        
        function self = set.password(self, newPassword)
            self.password = newPassword;
        end
        
        function self = set.server(self, newServer)
            self.server = newServer;
        end
        
        function self = set.recipients(self, newRecipients)
            self.recipients = newRecipients;
        end
        
        function self = set.subject(self, newSubject)
            self.subject = newSubject;
        end
        
        function self = set.message(self, newMessage)
            self.message = newMessage;
        end
    end
    
    methods
        function send(self)
            %
            % Send email message to recipients.
            %
            
            INTERNET_KEYWORD = 'Internet';
            E_MAIL_KEYWORD   = 'E_mail';
            SERVER_KEYWORD   = 'SMTP_Server';
            USERNAME_KEYWORD = 'SMTP_Username';
            PASSWORD_KEYWORD = 'SMTP_Password';
            AUTHENTICATION_KEYWORD      = 'mail.smtp.auth';
            TRUE_KEYWORD                = 'true';
            SOCKET_FACTORY_CLASS        = 'mail.smtp.socketFactory.class';
            SOCKET_JAVA_FACTORY_LIBRARY = 'javax.net.ssl.SSLSocketFactory';
            SOCKET_FACTORY_PORT         = 'mail.smtp.socketFactory.port';
            SOCKET_PORT_NUMBER          = '465';
            
            % Email can be sent if proper configurations have been made.
            if isempty(self.address) || isempty(self.password) || isempty(self.server)
                error = MException('email:send:noConfig', 'Email configurations is not performed...');
                error.throw();
            end
            
            % Email can be sent if there is a minimum of information
            % defined.
            if isempty(self.recipients) || isempty(self.subject)
                error = MException('email:send:noInfo', 'Recipients or subject not defined...');
                error.throw();
            end
            
            % Define settings
            setpref(INTERNET_KEYWORD,E_MAIL_KEYWORD  ,self.address );
            setpref(INTERNET_KEYWORD,SERVER_KEYWORD  ,self.server  );
            setpref(INTERNET_KEYWORD,USERNAME_KEYWORD,self.address );
            setpref(INTERNET_KEYWORD,PASSWORD_KEYWORD,self.password);

            props = java.lang.System.getProperties;
            props.setProperty(AUTHENTICATION_KEYWORD,TRUE_KEYWORD               );
            props.setProperty(SOCKET_FACTORY_CLASS  ,SOCKET_JAVA_FACTORY_LIBRARY);
            props.setProperty(SOCKET_FACTORY_PORT   ,SOCKET_PORT_NUMBER         );

            % Send email
            sendmail(self.recipients, self.subject, self.message);
        end
    end
    
    methods (Static)
        function self = EmailFromFile(filepath)
            %
            % Load email configuration from file information.
            %
            % INPUTS:
            %   1. filepath - local path of the file with the email
            %   configuration.
            %
            narginchk(1,1);
			
            % Create an empty email instance
            self = utils.Email();
            
            % Open file
            fileID   = fopen(filepath, 'r');
            
            % Load email configurations
            self.address  = fgetl(fileID);
            self.password = fgetl(fileID);
            self.server   = fgetl(fileID);
            
            % Close file
            fclose(fileID);
        end
    end
end