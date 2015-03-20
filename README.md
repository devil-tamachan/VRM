# VRM - Vector Rendering Method
A vector rendering script for Blender 2.73-

## Original version
[http://vrm.ao2.it/index.html](http://vrm.ao2.it/index.html) by Antonio Ospite

## libming (Support for Flash .swf)
[Precompiled Binary for Windows](https://github.com/devil-tamachan/libming/releases)

 - Blender2.73 32bit: Copy to C:\Program Files (x86)\Blender Foundation\Blender\2.73\python\lib
 - Blender2.73 64bit: Copy to C:\Program Files\Blender Foundation\Blender\2.73\python\lib

## TODO:
 - edge


# VRM
ベクトル化スクリプト for Blender 2.73～ 移植版

## オリジナルバージョン
[http://vrm.ao2.it/index.html](http://vrm.ao2.it/index.html) Antonio Ospite

## libming (Support for Flash .swf)
[Windows向けコンパイル済みバイナリ](https://github.com/devil-tamachan/libming/releases)

 - Blender2.73 32bit: 
  * 上記からWindows Python 3.4 32bit (Blender 2.73 32bit)ダウンロードして解凍。
  * 中身を C:\Program Files (x86)\Blender Foundation\Blender\2.73\python\lib へコピー
 - Blender2.73 64bit: 
  * 上記からWindows Python 3.4 64bit (Blender 2.73 64bit)ダウンロードして解凍。
  * 中身を C:\Program Files\Blender Foundation\Blender\2.73\python\lib へコピー

## 使用上のTips
 - 面は３角形に割られます
 - スムーズはかかりません（必要なら事前に適応してください）。あまり割りすぎると遅くなります。
 - ポリゴン数が多すぎると表示できないswfが出力される場合があります。
 - カメラ範囲からめちゃくちゃはみ出すデカポリゴンがあると稀にBlenderごと死ぬ場合があります。そういう場合は途中で割って小さくしてください
 - Flash Proに読み込み、毎フレームで"すべて選択"→"分解"すると劇的に軽量化できます。ただFlashProのバグ？で破損する場合があります。

## TODO:
 - edge出力未サポート