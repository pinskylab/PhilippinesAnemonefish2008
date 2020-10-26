These are the command lines to gather waypoint and track data from the 
GPS.  

gpsman -dev usb=/dev/ttyUSB0 getwrite WP GPX dump_waypoints.<date>.gpx
gpsman -dev usb=/dev/ttyUSB0 getwrite TR GPX dump_waypoints.<date>.gpx
