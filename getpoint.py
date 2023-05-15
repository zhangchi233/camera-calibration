import cv2


def mouse(event, x, y, flags, param):
    global flag, x1, y1, x2, y2, wx, wy, move_w, move_h, dst
    global zoom, zoom_w, zoom_h, img_zoom, flag_har, flag_var
    if event == cv2.EVENT_LBUTTONDOWN:
        xy = "%d,%d" % (x, y)
        a.append(x)
        b.append(y)
        cv2.circle(dst, (x, y), 1, (0, 0, 255), thickness=-1)
        cv2.putText(dst, xy, (x, y), cv2.FONT_HERSHEY_PLAIN,
                    1, (0, 0, 255), thickness=2)
        # cv2.imshow("image", img_original)
        print(x, y)
    # if event==cv2.EVENT_MOUSEWHEEL:

    cv2.imshow("img", dst)
    cv2.waitKey(1)


win_h, win_w = 600, 700  # 窗口宽高
wx, wy = 0, 0  # 窗口相对于原图的坐标
wheel_step, zoom = 0.05, 1  # 缩放系数， 缩放值
img_original = cv2.imread("ourcorner.png")  # 建议图片大于win_w*win_h(800*600)
dst = cv2.imread("ourcorner.png")
# dst = cv2.resize(img_original,None,fx=3,fy=3,interpolation=cv2.INTER_LINEAR)
# (img_original)
# cv2.imshow("imge",img_original)
a = []
b = []
cv2.namedWindow('img', cv2.WINDOW_NORMAL)
# cv2.moveWindow("img", 900, 800)

x, y = img_original.shape[0:2]

# dst = cv2.resize(img_original, (int(y * 2), int(x * 2)))
cv2.setMouseCallback('img', mouse)
cv2.waitKey()
cv2.destroyAllWindows()
#