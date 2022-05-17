# -*- coding: utf-8 -*-
#
#  Copyright 2022 Ramil Nugmanov <nougmanoff@protonmail.com>
#  This file is part of chython.
#
#  chython is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with this program; if not, see <https://www.gnu.org/licenses/>.
#
from asyncio import new_event_loop
from os.path import join
from tempfile import TemporaryDirectory


loop = browser = None


async def render(s, t, width, height, scale):
    page = await browser.newPage()
    await page.setViewport({'deviceScaleFactor': scale, 'width': width, 'height': height})
    await page.goto(f'file://{s}')
    element = await page.querySelector('svg')
    await element.screenshot({'path': t})
    await page.close()


def svg2png(svg: str, width: int = 1000, height: int = 1000, scale: float = 10.):
    global loop, browser

    if loop is None:  # lazy browser launcher
        from pyppeteer import launch

        loop = new_event_loop()
        browser = loop.run_until_complete(launch())
    elif browser is None:
        raise ImportError('pyppeteer initialization failed')

    with TemporaryDirectory() as tmpdir:
        with open(s := join(tmpdir, 'input.svg'), 'w') as f:
            f.write(svg)

        loop.run_until_complete(render(s, (t := join(tmpdir, 'output.png')), width, height, scale))

        with open(t, 'rb') as f:
            return f.read()


__all__ = ['svg2png']
