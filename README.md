php-astro
===========

Class Astro is small collection of the astronomical functions written in php. 

This class is based on the works:

Darin Koblick https://www.mathworks.com/matlabcentral/profile/authors/1284781-darin-koblick

Dubravko Loborec http://www.dubravkoloborec.com


## Installing

You can use [Composer](http://getcomposer.org/) to add the [package](https://packagist.org/packages/loborec/zf2-message) to your project:

```json
{
  "require": {
    "loborec/astro": "dev-master"
  }
}
```

## Functions
### Lunar Azimuth and Altitude Estimation Algorithm (LunarAzEl)
Predict the Lunar Azimuth and Altitude within +/- .2 deg of any lat and lon for a given UTC

This algorithm will accept a Latitude, Longitude and Altitude location as well as a specific universal coordinated time. It will use this information and calculate the position of the moon in a local coordinate frame (az and alt aka az and el).

Example Function Call: 
```php
$result = Astro::LunarAzEl(1991, 05, 19, 13, 0, 0, 50, 10, 0);
var_dump(result);
```
Result:
```txt
array (size=2)
  0 => float 110.63005216655
  1 => float 35.226177186756
```

Input List: 

UTC Date and Time - Year, Month, Day, Hour, Minute, Second

Latitude - Site Latitude in degrees -90:90 -> S(-) N(+) 

Longitude - Site Longitude in degrees -180:180 W(-) E(+) 

Altitude - Site Altitude in km

Output List: 

Az - Lunar Azimuth angle in degrees 

El - Lunar Elevation/Altitude Angle in degrees

Verified output by comparison with the following source data: http://aa.usno.navy.mil/data/docs/AltAz.php

Copyright (c) 2010, Darin Koblick

From MATLAB to php translated by Dubravko Loborec

### Julian date Algorithm (juliandate)
The Julian date (JD) of any instant is the Julian day number for the preceding noon in Universal Time plus the fraction of the day since that instant. Julian dates are expressed as a Julian day number with a decimal fraction added

Example converting date to UTC, and calculate juliandate;

```php
    $date = new DateTime('2017-07-20T10:20:00 CEST');
    $gmtime=$date->getTimestamp();

    $y=gmdate('Y', $gmtime);
    $m=gmdate('m', $gmtime);
    $d=gmdate('d', $gmtime);
    $h=gmdate('H', $gmtime);
    $mn=gmdate('i', $gmtime);
    $s=gmdate('s', $gmtime);

    echo juliandate($y, $m, $d, $h, $mn, $s);
```

Input: 

UTC Date and Time - Year, Month, Day, Hour, Minute, Second

Output:

Julian date