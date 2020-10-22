classdef kalman_filter < handle
   properties(Access = private)
      x
      x_dot
      P
      P_dot
      z
      H_K
      w
      Q
      K
   end
   methods
       function obj = kalman_filter()
           
       end
   end
   methods (Access = private)
       function project_state(obj)
           p = obj.x(1:3);
           b = obj.x(4:6);
           p_dot = 0.5 * (0.5 * (1 - p' * p) * eye(3) + tilde(p) + p * p') * (obj.w - b - eta_1);
           b_dot = eta_2;
           obj.x_dot = [p_dot;
                        b_dot];
           obj.x = obj.x + obj.x_dot * dt;
       end
       function project_error_covariance(obj)
           p = obj.x(1:3);
           df_dp = 0.5 * (p * obj.w' - obj.w * p' - tilde(obj.w) + (obj.w' * p) * eye(3));
           df_db = -0.5 * (0.5 * (1 - p'*p)*eye(3) + tilde(p) + p * p');
           F = [df_dp df_db;
               zeros(3) zeros(3)];
           G = [df_db zeros(3);
               zeros(3) eye(3)];
           obj.P_dot = F * obj.P + obj.P' * F + G * obj.Q * G';
           obj.P = obj.P + obj.P_dot * dt;
       end
       function kalman_gain(obj)
           p = obj.x(1:3);
           A = attitude_matrix();
           % TODO: Define BI ref. 1 MRPS
           %        Define R
           L = 4 / (1 + p' * p) ^ 2 * tilde (A * BI) * ((1 - p' * p) * eye(3) - 2 * tilde(p) + 2 * (p * p'));
           obj.H_K = [L zeros(3)];
           obj.K = (obj.P * obj.H_K') / (obj.H_K * obj.P * obj.H_K' + R);
       end
       function estimate_state(obj)
           % TODO: Define BI ref. 1 MRPS
           %        Define vk
           hk = attitude_matrix() * BI;
           zk = hk + vk;
           obj.x = obj.x + obj.K * (zk - hk);
       end
       function estimate_uncertainty(obj)
           obj.P = (eye(6) - obj.K * obj.H_K * obj.x) * obj.P;
       end
       function A = attitude_matrix(obj)
           p = obj.x(1:3);
           A = eye(3) - 4 * (1 - p' * p) / (1 + p' * p) ^ 2 * tilde(p) + 8 / (1 + p' * p) ^ 2 * tilde(p) ^ 2;
       end
   end
   methods (Static)
       function V = tilde(v)
           V = [0  -v(3) v(2);
               v(3) 0 -v(1);
               -v(2) v(1) 0];
       end
   end
end
