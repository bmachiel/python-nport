#!/bin/sh

hspice deemb_vandamme01.sp -o /tmp
hspice deemb_vandamme01_open.sp -o /tmp
hspice deemb_vandamme01_short1.sp -o /tmp
hspice deemb_vandamme01_short2.sp -o /tmp
hspice deemb_vandamme01_through.sp -o /tmp
