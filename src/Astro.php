<?php

    class Astro{

        /**
        * php translation of a cart2sph MATLAB function
        * [azimuth,elevation,r] = cart2sph(x,y,z)
        * 
        * @author Dubravko Loborec
        * @copyright  2017, Dubravko Loborec
        * @link http://www.dubravkoloborec.com
        * 
        * @param mixed $x
        * @param mixed $y
        * @param mixed $z
        * @return array[$azimuth, $elevation, $r]
        */
        public static function cart2sph($x,$y,$z){
            $azimuth = atan2($y,$x);
            $elevation = atan2($z, sqrt(pow($x,2) + pow($y,2)));
            $r = sqrt(pow($x,2) + pow($y,2) + pow($z,2));
            return [$azimuth, $elevation, $r];
        }

        /**
        * php translation of a sph2cart MATLAB function
        * [x,y,z] = sph2cart(azimuth,elevation,r)
        *  
        * @author Dubravko Loborec
        * @copyright  2017, Dubravko Loborec
        * @link http://www.dubravkoloborec.com
        * 
        * @param mixed $azimuth
        * @param mixed $elevation
        * @param mixed $r
        * @return array[$x, $y, $z]
        */
        public static function sph2cart($azimuth, $elevation, $r){
            $x = $r * cos($elevation) * cos($azimuth);
            $y = $r * cos($elevation) * sin($azimuth);
            $z = $r * sin($elevation);
            return [$x, $y, $z];
        }

        /**
        * php translation of a dot MATLAB function
        *      
        * @author Dubravko Loborec
        * @copyright  2017, Dubravko Loborec
        * @link http://www.dubravkoloborec.com
        * 
        * @param array $a
        * @param array $b
        * @return float
        */
        public static function dot($a, $b){
            $n=count($a);
            $sum=0;
            for ($x = 0; $x < $n; $x++) {
                $sum=$sum+$a[$x]*$b[$x];
            } 
            return $sum;   
        }

        /**
        * Returns integer value of the $value
        * 
        * @param mixed $value
        * @return int
        */
        public static function bigintval($value){
            $value = trim($value);
            if (ctype_digit($value)) {
                return $value;
            }
            $value = preg_replace("/[^0-9](.*)$/", '', $value);
            if (ctype_digit($value)) {
                return $value;
            }
            return 0;
        }

        /**
        * helper function
        * 
        * @param mixed $a
        * @param mixed $b
        * @param mixed $c
        */
        private static function sol_calc($a, $b, $c){
            $b[0]=$b[0]*$c;  
            $b[1]=$b[1]*$c;
            $b[2]=$b[2]*$c;  

            $x=$a[0][0]*$b[0]+$a[1][0]*$b[1]+$a[2][0]*$b[2];
            $y=$a[0][1]*$b[0]+$a[1][1]*$b[1]+$a[2][1]*$b[2];
            $z=$a[0][2]*$b[0]+$a[1][2]*$b[1]+$a[2][2]*$b[2];

            return [$x, $y, $z]; 
        }

        /**
        * Calculates juliandate from UTC date
        * php translation by Dubravko Loborec
        * 
        * @author Darin Koblick
        * @copyright  2010, Darin Koblick
        * @link https://www.mathworks.com/matlabcentral/fileexchange/22992-lunar-azimuth-and-altitude-estimation-algorithm
        * @license https://www.mathworks.com/matlabcentral/fileexchange/22992-lunar-azimuth-and-altitude-estimation-algorithm
        * 
        * @param integer $y //UTC year
        * @param integer $m //UTC month
        * @param integer $d //UTC year
        * @param integer $h //UTC hour
        * @param integer $mn //UTC minute
        * @param integer $s //UTC second
        * @return float
        */
        public static function juliandate($y, $m, $d, $h=0, $mn=0, $s=0){
            if ( $m <= 2 ){ // january & february
                $y  = $y - 1.0;
                $m =$m + 12.0;
            }          

            return self::bigintval( 365.25*($y + 4716.0)) + self::bigintval( 30.6001*( $m + 1.0)) + 2.0 - 
            self::bigintval( $y/100.0 ) + self::bigintval( self::bigintval( $y/100.0 )/4.0 ) + $d - 1524.5 + 
            ($h + $mn/60 + $s/3600)/24;
        }            

          /**
        * Predict the azimuth and elevation of the Sun within +/- 1 degree at any geodetic latitude, longitude and altitude.
        * 
        * php translation and modification by Dubravko Loborec
        * 
        * @author Darin Koblick
        * @copyright  2010, Darin Koblick
        * @link https://www.mathworks.com/matlabcentral/fileexchange/23051-vectorized-solar-azimuth-and-elevation-estimation
        * @license https://www.mathworks.com/matlabcentral/fileexchange/23051-vectorized-solar-azimuth-and-elevation-estimation
        * 
        * @param integer $y //UTC year
        * @param integer $m //UTC month
        * @param integer $d //UTC year
        * @param integer $h //UTC hour
        * @param integer $mn //UTC minute
        * @param integer $s //UTC second
        * @param float $Lat //Site Latitude in degrees -90:90 -> S(-) N(+) 
        * @param float $Lon //Site Longitude in degrees -180:180 W(-) E(+) 
        * @param float $Alt //Site Altitude in km
        * @return array[$Az, $h]
        */
        function SolarAzEl($y, $m, $d, $h=0, $mn=0, $s=0, $Lat, $Lon, $Alt){

            $jd=self::juliandate($y, $m, $d, $h, $mn, $s);

            $d = $jd-2451543.5;

            // Keplerian Elements for the Sun (geocentric)
            $w = 282.9404+4.70935e-5 * $d; //    (longitude of perihelion degrees)
            //a = 1.000000;//                  (mean distance, a.u.)
            $e = 0.016709-1.151e-9 * $d;//       (eccentricity)
            $M = fmod(356.0470+0.9856002585 * $d,360);//   (mean anomaly degrees)
            $L = $w + $M;                     //(Sun's mean longitude degrees)
            $oblecl = 23.4393-3.563e-7 * $d;  //(Sun's obliquity of the ecliptic)

            //auxiliary angle
            $EE = $M+(180/pi() ) * $e * sin($M * (pi() /180)) * (1+$e * cos($M * (pi() /180)));

            //rectangular coordinates in the plane of the ecliptic ($x axis toward
            //perhilion)
            $x = cos($EE * (pi() /180))-$e;
            $y = sin($EE * (pi() /180)) * sqrt(1-pow($e,2));

            //find the distance and true anomaly
            $r = sqrt(pow($x,2) + pow($y,2));
            $v = atan2($y,$x) * (180/pi() );

            //find the longitude of the sun
            $lon = $v + $w;

            //compute the ecliptic rectangular coordinates
            $xeclip = $r * cos($lon * (pi() /180));
            $yeclip = $r * sin($lon * (pi() /180));
            $zeclip = 0.0;

            //rotate these coordinates to equitorial rectangular coordinates
            $xequat = $xeclip;
            $yequat = $yeclip * cos($oblecl * (pi() /180))+$zeclip * sin($oblecl * (pi() /180));
            $zequat = $yeclip * sin(23.4406 * (pi() /180))+$zeclip * cos($oblecl * (pi() /180));

            //convert equatorial rectangular coordinates to $RA and Decl:
            $r = sqrt(pow($xequat,2) + pow($yequat,2) + pow($zequat,2))-($Alt/149598000); //roll up the altitude correction
            $RA = atan2($yequat,$xequat) * (180/pi() );
            $delta = asin($zequat/$r) * (180/pi() );

            //Following the $RA DEC to Az Alt conversion sequence explained here:
            //http://www.stargazing.net/kepler/altaz.html

            //Find the J2000 value
            //J2000 = $jd - 2451545.0;
            //hourvec = datevec(UTC);
            //UTH = hourvec(:,4) + hourvec(:,5)/60 + hourvec(:,6)/3600;
            $UTH=$h + $mn/60 + $s/3600;

            //Calculate local siderial time
            $GMST0=fmod($L+180,360)/15;
            $SIDTIME = $GMST0 + $UTH + $Lon/15;

            //Replace $RA with hour angle $HA
            $HA = ($SIDTIME * 15 - $RA);

            //convert to rectangular coordinate system
            $x = cos($HA * (pi() /180)) * cos($delta * (pi() /180));
            $y = sin($HA * (pi() /180)) * cos($delta * (pi() /180));
            $z = sin($delta * (pi() /180));

            //rotate this a$Long an axis going east-west.
            $xhor = $x * cos((90-$Lat) * (pi() /180))-$z * sin((90-$Lat) * (pi() /180));
            $yhor = $y;
            $zhor = $x * sin((90-$Lat) * (pi() /180))+$z * cos((90-$Lat) * (pi() /180));

            //Find the h and AZ 
            $Az = atan2($yhor,$xhor) * (180/pi() ) + 180;
            $El = asin($zhor) * (180/pi() );

            return [$Az, $El];

        }
        
        /**
        * Predict the Lunar Azimuth and Altitude within +/- .2 deg of any lat and lon for a given UTC
        * This algorithm will accept a Latitude, Longitude and Altitude location as well as a specific universal coordinated time. It will use this information and calculate the position of the moon in a local coordinate frame (az and alt aka az and el).
        * 
        * php translation and modification by Dubravko Loborec
        * 
        * @author Darin Koblick
        * @copyright  2010, Darin Koblick
        * @link https://www.mathworks.com/matlabcentral/fileexchange/22992-lunar-azimuth-and-altitude-estimation-algorithm
        * @license https://www.mathworks.com/matlabcentral/fileexchange/22992-lunar-azimuth-and-altitude-estimation-algorithm
        * 
        * @param integer $y //UTC year
        * @param integer $m //UTC month
        * @param integer $d //UTC year
        * @param integer $h //UTC hour
        * @param integer $mn //UTC minute
        * @param integer $s //UTC second
        * @param float $Lat //Site Latitude in degrees -90:90 -> S(-) N(+) 
        * @param float $Lon //Site Longitude in degrees -180:180 W(-) E(+) 
        * @param float $Alt //Site Altitude in km
        * @return array[$Az, $h]
        */
        public static function LunarAzEl($y, $m, $d, $h=0, $mn=0, $s=0, $Lat, $Lon, $Alt){

            while ($Lon > 180){
                $Lon = $Lon - 360;
            }
            while ($Lon < -180){
                $Lon = $Lon + 360; 
            }
            while ($Lat > 90){
                $Lat = $Lat - 360;
            }
            while ($Lat < -90){
                $Lat = $Lat + 360; 
            }

            $EarthRadEq = 6378.1370;

            $jd=self::juliandate($y, $m, $d, $h, $mn, $s);
            $d = $jd - 2451543.5;

            $N = 125.1228-0.0529538083*$d; //    (Long asc. node deg)
            $i = 5.1454; //                      (Inclination deg)
            $w = 318.0634 + 0.1643573223*$d; //  (Arg. of perigee deg)
            $a =  60.2666;//                     (Mean distance (Earth's Equitorial Radii)
            $e = 0.054900;//                     (Eccentricity)
            $M = fmod(115.3654+13.0649929509*$d,360);//    (Mean anomaly deg)


            $LMoon =  fmod($N + $w + $M,360);                 //(Moon's mean longitude deg)
            $FMoon =  fmod($LMoon - $N,360);                 //(Moon's argument of latitude)

            //Keplerian Elements of the Sun
            $wSun = fmod(282.9404 + 4.70935E-5*$d,360);    // (longitude of perihelion)
            $MSun = fmod(356.0470 + 0.9856002585*$d,360);  // (Sun mean anomaly)
            $LSun = fmod($wSun + $MSun,360);                 // (Sun's mean longitude)

            $DMoon =  $LMoon - $LSun;                     // (Moon's mean elongation)  


            //Calculate Lunar perturbations in Longitude
            $LunarPLon = [ 
                -1.274*sin(($M - 2 *$DMoon)*(pi() /180)), 
                .658*sin(2 *$DMoon*(pi() /180)), 
                -0.186*sin($MSun*(pi() /180)), 
                -0.059*sin((2 *$M-2 *$DMoon)*(pi() /180)), 
                -0.057*sin(($M-2*$DMoon + $MSun) *(pi() /180)), 
                .053 *sin(($M+2*$DMoon) *(pi() /180)), 
                .046 *sin((2*$DMoon-$MSun) *(pi() /180)), 
                .041 *sin(($M-$MSun) *(pi() /180)), 
                -0.035 *sin($DMoon *(pi() /180)),            
                -0.031 *sin(($M+$MSun) *(pi() /180)), 
                -0.015 *sin((2 *$FMoon-2 *$DMoon) *(pi() /180)), 
                .011 *sin(($M-4 *$DMoon) *(pi() /180))
            ];

            //Calculate Lunar perturbations in Latitude 
            $LunarPLat = [ 
                -0.173 *sin(($FMoon-2 * $DMoon) *(pi() /180)), 
                -0.055 *sin(($M-$FMoon-2 * $DMoon) *(pi() /180)), 
                -0.046 *sin(($M+$FMoon-2 * $DMoon) *(pi() /180)), 
                +0.033 *sin(($FMoon+2 * $DMoon) *(pi() /180)), 
                +0.017 *sin((2 * $M + $FMoon) *(pi() /180))
            ];

            //Calculate perturbations in Distance
            $LunarPDist = [                                   
                -0.58*cos(($M-2 * $DMoon) * (pi() /180)),
                -0.46 * cos(2 * $DMoon * (pi() /180))
            ];

            // Compute EE, the eccentric anomaly

            //$E0 is the eccentric anomaly approximation estimate 
            //(this will initially have $a relativly high error)
            $E0 = $M+(180./pi() ) * $e * sin($M * (pi() /180)) * (1+$e * cos($M * (pi() /180)));

            //Compute $E1 and set it to $E0 until the $E1 == $E0
            $E1 = $E0-($E0-(180/pi() ) * $e * sin($E0 * (pi() /180))-$M)/(1-$e*cos($E0 * (pi() /180)));

            while ($E1-$E0 > .000005){
                $E0 = $E1;
                $E1 = $E0-($E0-(180/pi() ) * $e * sin($E0 * (pi() /180))-$M)/(1-$e*cos($E0 * (pi() /180)));    
            }

            $EE = $E1;

            //Compute rectangular coordinates ($x,$y) in the plane of the lunar orbit
            $x = $a * (cos($EE * (pi() /180))-$e);
            $y = $a * sqrt(1-$e * $e) * sin($EE * (pi() /180));

            //convert this to distance and true anomaly
            $r = sqrt($x * $x + $y * $y);
            $v = atan2($y * (pi() /180),$x * (pi() /180)) * (180/pi() );

            //Compute moon's position in ecliptic coordinates
            $xeclip = $r * (cos($N * (pi() /180)) * cos(($v+$w) * (pi() /180))-sin($N * (pi() /180)) * sin(($v+$w) * (pi() /180)) * cos($i * (pi() /180)));
            $yeclip = $r * (sin($N * (pi() /180)) * cos(($v+$w) * (pi() /180))+cos($N * (pi() /180))*sin((($v+$w) * (pi() /180)))*cos($i * (pi() /180)));
            $zeclip = $r * sin(($v+$w) * (pi() /180)) * sin($i * (pi() /180));

            //Add the calculated lunar perturbation terms to increase fmodel fidelity
            list($eLon, $eLat, $eDist) = self::cart2sph($xeclip,$yeclip,$zeclip);

            list($xeclip, $yeclip, $zeclip) = self::sph2cart($eLon + array_sum($LunarPLon) * (pi() /180), 
                $eLat + array_sum($LunarPLat) * (pi() /180), 
                $eDist + array_sum($LunarPDist));


            //  clear $eLon $eLat $eDist;
            $eLon=0;
            $eLat=0;
            $eDist=0;

            //convert the latitude and longitude to right ascension RA and declination
            //delta
            $T = ($jd-2451545.0)/36525.0;

            //Generate $a rotation matrix for ecliptic to equitorial
            //RotM=rotm_coo('$e',$jd);
            //See rotm_coo.$M for obl and rotational matrix transformation
            $Obl = 23.439291 - 0.0130042 * $T - 0.00000016 * $T * $T + 0.000000504 * $T * $T * $T;
            $Obl = $Obl * (pi() /180);

            //$RotM = [1 0 0; 0 cos($Obl) sin($Obl); 0 -sin($Obl) cos($Obl)]'; /////////////////////

            $RotM= [[1, 0, 0],[ 0, cos($Obl), sin($Obl)],[0, -sin($Obl), cos($Obl)]];


            //Apply the rotational matrix to the ecliptic rectangular coordinates
            //Also, convert units to km instead of earth equatorial radii   

            $sol=self::sol_calc($RotM, [$xeclip, $yeclip, $zeclip], $EarthRadEq);

            //Find the equatorial rectangular coordinates of the location specified
            list($xel, $yel, $zel) = self::sph2cart($Lon * (pi() /180),$Lat * (pi() /180),$Alt+$EarthRadEq);

            //Find the equatorial rectangular coordinates of the location @ sea level
            list($xsl, $ysl, $zsl) = self::sph2cart($Lon * (pi() /180),$Lat * (pi() /180),$EarthRadEq);

            //Find the Angle Between sea level coordinate vector and the moon vector
            /*$theta1 = 180 - acosd(dot([$xsl $ysl $zsl], [$sol(1)-$xsl $sol(2)-$ysl $sol(3)-$zsl]) 
            /(sqrt($xsl.^2 + $ysl.^2 + $zsl.^2) 
            * sqrt(($sol(1)-$xsl).^2 + ($sol(2)-$ysl).^2 + ($sol(3)-$zsl).^2)));*/

            $theta1 = 180 - rad2deg(acos(
                self::dot( [$xsl, $ysl, $zsl], [$sol[0]-$xsl, $sol[1]-$ysl, $sol[2]-$zsl] )
                /(sqrt(pow($xsl,2) + pow($ysl,2) + pow($zsl,2)) 
                    * sqrt( pow($sol[0]-$xsl,2) + pow($sol[1]-$ysl,2) + pow($sol[2]-$zsl,2)))

            ));         

            //Find the Angle Between the same coordinates but at the specified elevation
            /*$theta2 = 180 - acosd(dot([$xel $yel $zel],[$sol(1)-$xel $sol(2)-$yel $sol(3)-$zel]) ...
            ./(sqrt($xel.^2 + $yel.^2 + $zel.^2) ...
            * sqrt(($sol(1)-$xel).^2 + ($sol(2)-$yel).^2 + ($sol(3)-$zel).^2)));*/

            $theta2 = 180 - rad2deg(acos(
                self::dot( [$xel, $yel, $zel], [$sol[0]-$xel, $sol[1]-$yel, $sol[2]-$zel] )
                /(sqrt(pow($xel,2) + pow($yel,2) + pow($zel,2)) 
                    *sqrt( pow($sol[0]-$xel,2) + pow($sol[1]-$yel,2) + pow($sol[2]-$zel,2)))
            ));     

            //Find the Difference Between the two angles (+|-) is important
            $thetaDiff = $theta2 - $theta1;

            // equatorial to horizon coordinate t$RAnsformation
            list($RA, $delta) = self::cart2sph($sol[0],$sol[1],$sol[2]);
            $delta = $delta * (180/pi() );
            $RA = $RA * (180/pi() );

            //Following the $RA DEC to Az $Alt conversion sequence explained here:
            //http://www.stargazing.net/kepler/altaz.html

            //Find the $J2000 value
            $J2000 = $jd - 2451545.0;

            //[Y,M,D,H,MN,S]
            //$hourvec = datevec($UTC, 'yyyy/mm/dd HH:MM:SS');
            //$UTH = hourvec(4) + hourvec(5)/60 + hourvec(6)/3600;

            $UTH=$h + $mn/60 + $s/3600;

            var_dump($mn/60);
            //Calculate local siderial time
            $LST = fmod(100.46+0.985647 * $J2000+$Lon+15*$UTH,360);

            //Replace $RA with hour angle HA
            $HA = $LST-$RA;

            //Find t$he $h and $Az at t$he current $LST
            $h = asin(sin($delta * (pi() /180)) * sin($Lat * (pi() /180)) + cos($delta * (pi() /180)) * cos($Lat * (pi() /180)) * cos($HA * (pi() /180))) * (180/pi() );
            $Az = acos((sin($delta * (pi() /180)) - sin($h * (pi() /180)) * sin($Lat * (pi() /180)))/(cos($h * (pi() /180)) * cos($Lat * (pi() /180)))) * (180/pi() );

            //Add in t$he angle offset due to t$he specified site elevation
            $h = $h + $thetaDiff;

            if (sin($HA * (pi() /180)) >= 0){
                $Az = 360-$Az; 
            }

            //Apply Pa$RAlax Correction if we are still on eart$h
            if ($Alt < 100){
                $horParal = 8.794/($r*6379.14/149.59787e6);
                $p = asin(cos($h * (pi() /180))*sin(($horParal/3600) * (pi() /180))) * (180/pi() );
                $h = $h-$p;
            }
            return [$Az, $h];
        }

      

    }

