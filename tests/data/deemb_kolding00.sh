#!/bin/sh

hspice deemb_kolding00.sp -o /tmp
hspice deemb_kolding00_simple_open.sp -o /tmp
hspice deemb_kolding00_simple_short.sp -o /tmp
hspice deemb_kolding00_open.sp -o /tmp
hspice deemb_kolding00_short1.sp -o /tmp
hspice deemb_kolding00_short2.sp -o /tmp
