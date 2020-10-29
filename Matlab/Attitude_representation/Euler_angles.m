classdef Euler_angles
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods (Static)
        function DCM = getDCM(angles, group, rad)
            if (nargin<3), or isempty(rad)
              rad = false;
            end
            if (~ rad)
                angles = angles / (180 / pi);
            end
            if group == "121"
                DCM = c121(angles);
            elseif group == "123"
                DCM = c123(angles);
            elseif group == "131"
                DCM = c131(angles);
            elseif group == "132"
                DCM = c132(angles);
            elseif group == "212"
                DCM = c212(angles);
            elseif group == "213"
                DCM = c213(angles);
            elseif group == "231"
                DCM = c231(angles);
            elseif group == "232"
                DCM = c232(angles);
            elseif group == "312"
                DCM = c312(angles);
            elseif group == "313"
                DCM = c313(angles);
            elseif group == "321"
                DCM = c321(angles);
            elseif group == "323"
                DCM = c323(angles);
            else
                disp("Error. There is no", group, "group")
                DCM = -1;
            end
        end
        function angles = get_Euler_angles(DCM, group, rad)
            if (nargin<3), or isempty(rad)
              rad = false;
            end
            if group == "121"
                angles = angles121(DCM);
            elseif group == "123"
                angles = angles123(DCM);
            elseif group == "131"
                angles = angles131(DCM);
            elseif group == "132"
                angles = angles132(DCM);
            elseif group == "212"
                angles = angles212(DCM);
            elseif group == "213"
                angles = angles213(DCM);
            elseif group == "231"
                angles = angles231(DCM);
            elseif group == "232"
                angles = angles232(DCM);
            elseif group == "312"
                angles = angles312(DCM);
            elseif group == "313"
                angles = angles313(DCM);
            elseif group == "321"
                angles = angles321(DCM);
            elseif group == "323"
                angles = angles323(DCM);
            else
                disp("Error. There is no", group, "group")
                angles = -1;
            end
            if (~ rad)
                angles = angles * (180 / pi);
            end
        end
        function rates = get_rates(angles, w, group, rad)
            if (nargin<3), or isempty(rad)
              rad = false;
            end
            if (~ rad)
                angles = angles / (180 / pi);
            end
            if group == "121"
                rates = rates121(angles) * w;
            elseif group == "123"
                rates = rates123(angles) * w;
            elseif group == "131"
                rates = rates131(angles) * w;
            elseif group == "132"
                rates = rates132(angles) * w;
            elseif group == "212"
                rates = rates212(angles) * w;
            elseif group == "213"
                rates = rates213(angles) * w;
            elseif group == "231"
                rates = rates231(angles) * w;
            elseif group == "232"
                rates = rates232(angles) * w;
            elseif group == "312"
                rates = rates312(angles) * w;
            elseif group == "313"
                rates = rates313(angles) * w;
            elseif group == "321"
                rates = rates321(angles) * w;
            elseif group == "323"
                rates = rates323(angles) * w;
            else
                disp("Error. There is no", group, "group")
                rates = -1;
            end
        end
        function angles = add_rotations(rot1, rot2, group1, group2, rad)
            if (nargin<3), or isempty(rad)
              rad = false;
            end
            if (~ rad)
                rot1 = rot1 / (180 / pi);
                rot2 = rot2 / (180 / pi);
            end
            dcm1 = getDCM(rot1, group1);
            dcm2 = getDCM(rot2, group2);
            dcm = DCM.addition(dcm1, dcm2);
            angles = get_euler_angles(dcm, group1, rad);
        end
        function dcm = c121(angles)
            dcm = [
                cos(angles(2)),
            sin(angles(2)) * sin(angles(1)),
            -sin(angles(2)) * cos(angles(1));
            sin(angles(3)) * sin(angles(2)),
            -sin(angles(3)) * cos(angles(2)) * sin(angles(1)) + cos(angles(3)) * cos(angles(1)),
            sin(angles(3)) * cos(angles(2)) * cos(angles(1)) + cos(angles(3)) * sin(angles(1));
            cos(angles(3)) * sin(angles(2)),
            -cos(angles(3)) * cos(angles(2)) * sin(angles(1)) - sin(angles(3)) * cos(angles(1)),
            cos(angles(3)) * cos(angles(2)) * cos(angles(1)) - sin(angles(3)) * sin(angles(1))
            ];
        end
        function dcm = c123(angles)
            dcm = [
                cos(angles(3)) * cos(angles(2)),
                cos(angles(3)) * sin(angles(2)) * sin(angles(1)) + sin(angles(3)) * cos(angles(1)),
                -cos(angles(3)) * sin(angles(2)) * cos(angles(1)) + sin(angles(3)) * sin(angles(1));
                -sin(angles(3)) * cos(angles(2)),
                -sin(angles(3)) * sin(angles(2)) * sin(angles(1)) + cos(angles(3)) * cos(angles(1)),
                sin(angles(3)) * sin(angles(2)) * cos(angles(1)) + cos(angles(3)) * sin(angles(1));
                sin(angles(2)),
                -cos(angles(2)) * sin(angles(1)),
                cos(angles(2)) * cos(angles(1))
            ];
        end
        function dcm = c131(angles)
            dcm = [
                cos(angles(2)),
                sin(angles(2)) * cos(angles(1)),
                sin(angles(2)) * sin(angles(1));
                -cos(angles(3)) * sin(angles(2)),
                cos(angles(3)) * cos(angles(2)) * cos(angles(1)) - sin(angles(3)) * sin(angles(1)),
                cos(angles(3)) * cos(angles(2)) * sin(angles(1)) + sin(angles(3)) * cos(angles(1));
                sin(angles(3)) * sin(angles(2)),
                -sin(angles(3)) * cos(angles(2)) * cos(angles(1)) - cos(angles(3)) * sin(angles(1)),
                -sin(angles(3)) * cos(angles(2)) * sin(angles(1)) + cos(angles(3)) * cos(angles(1))
            ];
        end
        function dcm = c132(angles)
            dcm = [
                cos(angles(3)) * cos(angles(2)),
                cos(angles(3)) * sin(angles(2)) * cos(angles(1)) + sin(angles(3)) * sin(angles(1)),
                cos(angles(3)) * sin(angles(2)) * sin(angles(1)) - sin(angles(3)) * cos(angles(1));
                -sin(angles(2)),
                cos(angles(2)) * cos(angles(1)),
                cos(angles(2)) * sin(angles(1));
                sin(angles(3)) * cos(angles(2)),
                sin(angles(3)) * sin(angles(2)) * cos(angles(1)) - cos(angles(3)) * sin(angles(1)),
                sin(angles(3)) * sin(angles(2)) * sin(angles(1)) + cos(angles(3)) * cos(angles(1))
            ];
        end
        function dcm = c212(angles)
            dcm = [
                -sin(angles(3)) * cos(angles(2)) * sin(angles(1)) + cos(angles(3)) * cos(angles(1)),
                sin(angles(3)) * sin(angles(2)),
                -sin(angles(3)) * cos(angles(2)) * cos(angles(1)) - cos(angles(3)) * sin(angles(1));
                sin(angles(2)) * sin(angles(1)),
                cos(angles(2)),
                sin(angles(2)) * cos(angles(1));
                cos(angles(3)) * cos(angles(2)) * sin(angles(1)) + sin(angles(3)) * cos(angles(1)),
                -cos(angles(3)) * sin(angles(2)),
                cos(angles(3)) * cos(angles(2)) * cos(angles(1)) - sin(angles(3)) * sin(angles(1));
            ];
        end
        function dcm = c213(angles)
            dcm = [
                sin(angles(3)) * sin(angles(2)) * sin(angles(1)) + cos(angles(3)) * cos(angles(1)),
                sin(angles(3)) * cos(angles(2)),
                sin(angles(3)) * sin(angles(2)) * cos(angles(1)) - cos(angles(3)) * sin(angles(1));
                cos(angles(3)) * sin(angles(2)) * sin(angles(1)) - sin(angles(3)) * cos(angles(1)),
                cos(angles(3)) * cos(angles(2)),
                cos(angles(3)) * sin(angles(2)) * cos(angles(1)) + sin(angles(3)) * sin(angles(1));
                cos(angles(2)) * cos(angles(1)),
                -sin(angles(2)),
                cos(angles(2)) * cos(angles(1));
            ];
        end
        function dcm = c231(angles)
            dcm = [
                cos(angles(2)) * cos(angles(1)),
                sin(angles(2)),
                -cos(angles(2)) * sin(angles(1));
                -cos(angles(3)) * sin(angles(2)) * cos(angles(1)) + sin(angles(3)) * sin(angles(1)),
                cos(angles(3)) * cos(angles(2)),
                cos(angles(3)) * sin(angles(2)) * sin(angles(1)) + sin(angles(3)) * cos(angles(1));
                sin(angles(3)) * sin(angles(2)) * cos(angles(1)) + cos(angles(3)) * sin(angles(1)),
                -sin(angles(3)) * cos(angles(2)),
                -sin(angles(3)) * sin(angles(2)) * sin(angles(1)) + cos(angles(3)) * cos(angles(1))
            ];
        end
        function dcm = c232(angles)
            dcm = [
                cos(angles(3)) * cos(angles(2)) * cos(angles(1)) - sin(angles(3)) * sin(angles(1)),
                cos(angles(3)) * sin(angles(2)),
                -cos(angles(3)) * cos(angles(2)) * sin(angles(1)) - sin(angles(3)) * cos(angles(1));
                -sin(angles(2)) * cos(angles(1)),
                cos(angles(2)),
                sin(angles(2)) * sin(angles(1));
                sin(angles(3)) * cos(angles(2)) * cos(angles(1)) + cos(angles(3)) * sin(angles(1)),
                sin(angles(3)) * sin(angles(2)),
                -sin(angles(3)) * cos(angles(2)) * sin(angles(1)) + cos(angles(3)) * cos(angles(1))
            ];
        end
        function dcm = c312(angles)
            dcm = [
                -sin(angles(3)) * sin(angles(2)) * sin(angles(1)) + cos(angles(3)) * cos(angles(1)),
                sin(angles(3)) * sin(angles(2)) * cos(angles(1)) + cos(angles(3)) * sin(angles(1)),
                -sin(angles(3)) * cos(angles(2));
                -cos(angles(2)) * sin(angles(1)),
                cos(angles(2)) * cos(angles(1)),
                sin(angles(2));
                cos(angles(3)) * sin(angles(2)) * sin(angles(1)) + sin(angles(3)) * cos(angles(1)),
                -cos(angles(3)) * sin(angles(2)) * cos(angles(1)) + sin(angles(3)) * sin(angles(1)),
                cos(angles(3)) * cos(angles(2))
            ];
        end
        function dcm = c313(angles)
            dcm = [
                -sin(angles(3)) * cos(angles(2)) * sin(angles(1)) + cos(angles(3)) * cos(angles(1)),
                sin(angles(3)) * cos(angles(2)) * cos(angles(1)) + cos(angles(3)) * sin(angles(1)),
                sin(angles(3)) * sin(angles(2));
                -cos(angles(3)) * cos(angles(2)) * sin(angles(1)) - sin(angles(3)) * cos(angles(1)),
                cos(angles(3)) * cos(angles(2)) * cos(angles(1)) - sin(angles(3)) * sin(angles(1)),
                cos(angles(3)) * sin(angles(2));
                sin(angles(2)) * sin(angles(1)),
                -sin(angles(2)) * cos(angles(1)),
                cos(angles(2))
            ];
        end
        function dcm = c321(angles)
            dcm = [
                cos(angles(2)) * cos(angles(1)),
                cos(angles(2)) * sin(angles(1)),
                -sin(angles(2));
                sin(angles(3)) * sin(angles(2)) * cos(angles(1)) - cos(angles(3)) * sin(angles(1)),
                sin(angles(3)) * sin(angles(2)) * sin(angles(1)) + cos(angles(3)) * cos(angles(1)),
                sin(angles(3)) * cos(angles(2));
                cos(angles(3)) * sin(angles(2)) * cos(angles(1)) + sin(angles(3)) * sin(angles(1)),
                cos(angles(3)) * sin(angles(2)) * sin(angles(1)) - sin(angles(3)) * cos(angles(1)),
                cos(angles(3)) * cos(angles(2))
            ];
        end
        function dcm = c323(angles)
            dcm = [
                cos(angles(3)) * cos(angles(2)) * cos(angles(1)) - sin(angles(3)) * sin(angles(1)),
                cos(angles(3)) * cos(angles(2)) * sin(angles(1)) + sin(angles(3)) * cos(angles(1)),
                -cos(angles(3)) * sin(angles(2));
                -sin(angles(3)) * cos(angles(2)) * cos(angles(1)) - cos(angles(3)) * sin(angles(1)),
                -sin(angles(3)) * cos(angles(2)) * sin(angles(1)) + cos(angles(3)) * cos(angles(1)),
                sin(angles(3)) * sin(angles(2));
                sin(angles(2)) * cos(angles(1)),
                sin(angles(2)) * sin(angles(1)),
                cos(angles(2))
            ];
        end
        function angles = angles121(dcm)
            angles = [
                arctan(dcm(1, 2), -dcm(1, 3));
                acos(dcm(1, 1));
                arctan(dcm(2, 1), dcm(3, 1))
            ];
        end
        function angles = angles123(dcm)
            angles = [
                arctan(-dcm(3, 2), dcm(3, 3));
                asin(dcm(3, 1));
                arctan(-dcm(2, 1), dcm(1, 1))
            ];
        end
        function angles = angles131(dcm)
            angles = [
                arctan(dcm(1, 3), dcm(1, 2));
                acos(dcm(1, 1));
                arctan(dcm(3, 1), -dcm(2, 1))
            ];
        end
        function angles = angles132(dcm)
            angles = [
                arctan(dcm(2, 3), dcm(2, 2));
                asin(-dcm(1, 1));
                arctan(dcm(3, 1), dcm(1, 1))
            ];
        end
        function angles = angles212(dcm)
            angles = [
                arctan(dcm(2, 1), dcm(2, 3));
                acos(dcm(2, 2));
                arctan(dcm(1, 2), -dcm(3, 2))
            ];
        end
        function angles = angles213(dcm)
            angles = [
                arctan(dcm(3, 1), dcm(3, 3));
                -asin(dcm(3, 2));
                arctan(dcm(1, 2), dcm(2, 2))
            ];
        end
        function angles = angles231(dcm)
            angles = [
                arctan(-dcm(1, 3), dcm(1, 1));
                asin(dcm(1, 2));
                arctan(-dcm(3, 2), dcm(2, 2))
            ];
        end
        function angles = angles232(dcm)
            angles = [
                arctan(dcm(2, 3), -dcm(2, 1));
                acos(dcm(2, 2));
                arctan(dcm(3, 2), dcm(1, 2))
            ];
        end
        function angles = angles312(dcm)
            angles = [
                arctan(dcm(2, 2), -dcm(2, 1));
                asin(dcm(2, 3));
                arctan(-dcm(1, 3), dcm(3, 3))
            ];
        end
        function angles = angles313(dcm)
            angles = [
                arctan(dcm(3, 1), -dcm(3, 2));
                acos(dcm(3, 3));
                arctan(dcm(1, 3), dcm(2, 3))
            ];
        end
        function angles = angles321(dcm)
            angles = [
                arctan(dcm(1, 2), dcm(1, 1));
                -asin(dcm(1, 3));
                arctan(dcm(2, 3), dcm(3, 3))
            ];
        end
        function angles = angles323(dcm)
            angles = [
                arctan(dcm(3, 2), dcm(3, 1));
                acos(dcm(3, 3));
                arctan(dcm(2, 3), -dcm(1, 3))
            ];
        end
        function rates = rates121(angles)
            rates = [
                0, sin(angles(3)), cos(angles(3));
                0, sin(angles(2)) * cos(angles(3)), - sin(angles(2)) * sin(angles(3));
                sin(angles(2)), -cos(angles(2)) * sin(angles(3)), -cos(angles(2)) * cos(angles(3))
            ]/sin(angles(2));
        end
        function rates = rates123(angles)
            rates = [
                cos(angles(3)), -sin(angles(3)), 0;
                cos(angles(2)) * sin(angles(3)), cos(angles(2)) * cos(angles(3)), 0;
                -sin(angles(2)) * cos(angles(3)), sin(angles(2)) * sin(angles(3)), cos(angles(2))
            ]/cos(angles(2));
        end
        function rates = rates131(angles)
            rates = [
                0, -cos(angles(3)), sin(angles(3));
                0, sin(angles(2)) * sin(angles(3)), sin(angles(2)) * cos(angles(3));
                sin(angles(2)), cos(angles(2)) * cos(angles(3)), - cos(angles(2)) * sin(angles(3))
            ]/sin(angles(2));
        end
        function rates = rates132(angles)
            rates = [
                cos(angles(3)), 0, sin(angles(3));
                -cos(angles(2)) * sin(angles(3)), 0, cos(angles(2)) * cos(angles(3));
                sin(angles(2)) * cos(angles(3)), cos(angles(2)), sin(angles(2)) * sin(angles(3))
            ]/cos(angles(2));
        end
        function rates = rates212(angles)
            rates = [
                sin(angles(3)), 0 - cos(angles(3));
                sin(angles(2)) * cos(angles(3)), 0, sin(angles(2)) * sin(angles(3));
                -cos(angles(2)) * sin(angles(3)), sin(angles(2)), cos(angles(2)) * cos(angles(3))
            ]/sin(angles(2));
        end
        function rates = rates213(angles)
            rates = [
                sin(angles(3)), cos(angles(3)), 0;
                cos(angles(2)) * cos(angles(3)), -cos(angles(2)) * sin(angles(3)), 0;
                sin(angles(2)) * sin(angles(3)), sin(angles(2)) * cos(angles(3)), cos(angles(2))
            ]/cos(angles(2));
        end
        function rates = rates231(angles)
            rates = [
                0, cos(angles(3)), -sin(angles(3));
                0, cos(angles(2)) * sin(angles(3)), cos(angles(2)) * cos(angles(3));
                cos(angles(2)), -sin(angles(2)) * cos(angles(3)), sin(angles(2)) * sin(angles(3))
            ]/cos(angles(2));
        end
        function rates = rates232(angles)
            rates = [
                cos(angles(3)), 0, sin(angles(3));
                -sin(angles(2)) * sin(angles(3)), 0, sin(angles(2)) * cos(angles(3));
                -cos(angles(2)) * cos(angles(3)), sin(angles(2)), -cos(angles(2)) * sin(angles(3))
            ]/sin(angles(2));
        end
        function rates = rates312(angles)
            rates = [
                -sin(angles(3)), 0, cos(angles(3));
                cos(angles(2)) * cos(angles(3)), 0, cos(angles(2)) * sin(angles(3));
                sin(angles(2)) * sin(angles(3)), cos(angles(2)), -sin(angles(2)) * cos(angles(3))
            ]/cos(angles(2));
        end
        function rates = rates313(angles)
            rates = [
                sin(angles(3)), cos(angles(3)), 0;
                sin(angles(2)) * cos(angles(3)), -sin(angles(2)) * sin(angles(3)), 0;
                -cos(angles(2)) * sin(angles(3)), -cos(angles(2)) * cos(angles(3)), sin(angles(2))
            ]/sin(angles(2));
        end
        function rates = rates321(angles)
            rates = [
                0.0, sin(angles(3)), cos(angles(3));
                0.0, cos(angles(2)) * cos(angles(3)), -cos(angles(2)) * sin(angles(3));
                cos(angles(2)), sin(angles(2)) * sin(angles(3)), sin(angles(2)) * cos(angles(3))
            ]/cos(angles(2));
        end
        function rates = rates323(angles)
            rates = [
                -cos(angles(3)), sin(angles(3)), 0;
                sin(angles(2)) * sin(angles(3)), sin(angles(2)) * cos(angles(3)), 0;
                cos(angles(2)) * cos(angles(3)), -cos(angles(2)) * sin(angles(3)), sin(angles(2))
            ]/sin(angles(2));
        end
        
        
    end
end

